"""
Step: Gaia Validation.

For each tile:
  1. Loads the SExtractor catalog and selects PSF-quality sources
     (CLASS_STAR >= 0.84, FLAGS <= 7).
  2. Splits off Outer_2 (the faintest region) as SExtractor-only —
     these stars are too faint/numerous to rely on Gaia for.
  3. Queries Gaia DR3 within the tile's sky footprint (cached to disk
     per tile so repeat runs don't re-query), restricted to sources
     with G >= 11 and either a high stellar classification probability
     or G < 14.
  4. Cross-matches SExtractor sources to Gaia via a kd-tree nearest
     neighbor search on the unit sphere, with a per-tile match radius
     found from the separation histogram (not a fixed arcsec cut).
  5. Combines the Gaia-matched and SE-only (Outer_2) sources, drops
     duplicates, and rejects sources whose brightest pixel isn't close
     to the cutout center (saturated sources are exempted from this
     check, since their peak can be genuinely off-center from bleed).
  6. Labels every surviving star's REGION via magnitude alone
     (config.MAG_CUTS) — this is the authoritative per-star patch
     assignment; see the note in config.py about reconciling this
     with cutouts.py's flux-radius-based PATCHES.

Saves the combined master.csv plus one matched CSV per tile, all
under config.GAIA_OUTPUT_DIR. The raw per-tile Gaia query results
are cached separately under config.GAIA_CACHE_DIR (output/gaia_cache/),
auto-fetched on first run — there's nothing to provide as input.
"""

import math
import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
from astroquery.gaia import Gaia
from scipy.spatial import cKDTree

from . import config
from .helpers import finalize_plot, setup_plot_mode
from .point_sources import _read_sextractor_catalog


def _xyz(ra, dec):
    """Convert RA/Dec (degrees) to unit-sphere Cartesian coordinates, for kd-tree matching."""
    ra, dec = np.deg2rad(ra), np.deg2rad(dec)
    return np.column_stack([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)])


def _is_centered(data, x, y, box=21, tol=3):
    """
    True if the brightest pixel within a `box`x`box` cutout at (x, y)
    is within `tol` pixels of the cutout center. Returns False (not
    True) if the requested box falls off the edge of `data`.
    """
    h = box // 2
    xi, yi = int(round(x)), int(round(y))
    cut = data[yi - h:yi + h + 1, xi - h:xi + h + 1]
    if cut.shape != (box, box):
        return False
    py, px = np.unravel_index(np.nanargmax(cut), cut.shape)
    return np.sqrt((px - h) ** 2 + (py - h) ** 2) <= tol


def get_region(gmag, mag_cal, mag_cuts=None):
    """
    Region label for a star, using Gaia G magnitude when available and
    falling back to the calibrated SExtractor magnitude otherwise.
    """
    mag_cuts = mag_cuts or config.MAG_CUTS
    mag = gmag if not np.isnan(gmag) else mag_cal
    for region, (lo, hi) in mag_cuts.items():
        if lo <= mag < hi:
            return region
    return None


def cross_match_tile(tile_name, mag_cuts=None, zp=None, use_cache=True):
    """
    Run the full Gaia cross-match + centering-check pipeline for one
    tile. Returns (fits_path, matched_df); does not write to disk
    itself — build_master_catalog() saves the matched CSV under
    config.GAIA_OUTPUT_DIR once all tiles are combined.
    """
    mag_cuts = mag_cuts or config.MAG_CUTS
    zp = config.GAIA_ZP if zp is None else zp

    print(f"\n{'=' * 50}  {tile_name}")

    fits_path = os.path.join(config.FITS_DIR, f"{tile_name}.fits")
    cat_file = os.path.join(config.SE_IO_DIR, f"outparams_{tile_name}.cat")
    gaia_dir = config.tile_gaia_dir(tile_name)
    os.makedirs(gaia_dir, exist_ok=True)

    with fits.open(fits_path, memmap=True) as hdul:
        wcs = WCS(hdul[0].header)
        header = hdul[0].header
        nx, ny = header["NAXIS1"], header["NAXIS2"]

        catalog = _read_sextractor_catalog(cat_file)
        sky = wcs.pixel_to_world(catalog["X_IMAGE"].values, catalog["Y_IMAGE"].values)
        catalog["RA"] = sky.ra.deg
        catalog["Dec"] = sky.dec.deg
        catalog["MAG_CAL"] = catalog["MAG_AUTO"] + zp

        catalog_psf = catalog[
            (catalog["CLASS_STAR"] >= 0.84) & (catalog["FLAGS"] <= 7)
        ].copy().reset_index(drop=True)
        print(f"  PSF sources (CLASS_STAR\u22650.84, FLAGS\u22647): {len(catalog_psf)}")

        # Outer_2 — SExtractor only, no Gaia cross-match attempted
        o2_lo, o2_hi = mag_cuts["Outer_2"]
        se_outer2 = catalog_psf[
            (catalog_psf["MAG_CAL"] >= o2_lo) & (catalog_psf["MAG_CAL"] < o2_hi)
        ].copy().reset_index(drop=True)
        for col in ["GAIA_RA", "GAIA_Dec", "GAIA_Gmag", "Sep_arcsec"]:
            se_outer2[col] = np.nan
        se_outer2["TILE"] = tile_name
        se_outer2["SOURCE"] = "SE_only"
        print(f"  Outer_2 SE-only: {len(se_outer2)}")

        # Gaia query (cached per tile)
        gaia_out = os.path.join(gaia_dir, "gaia_sources.csv")
        corners_sky = [wcs.pixel_to_world(x, y) for x, y in [(0, 0), (nx, 0), (0, ny), (nx, ny)]]
        ra_min, ra_max = min(c.ra.deg for c in corners_sky), max(c.ra.deg for c in corners_sky)
        dec_min, dec_max = min(c.dec.deg for c in corners_sky), max(c.dec.deg for c in corners_sky)

        if use_cache and os.path.exists(gaia_out):
            gaia = pd.read_csv(gaia_out).astype({"ra": float, "dec": float})
            print(f"  Gaia (cached): {len(gaia)}")
        else:
            Gaia.ROW_LIMIT = -1
            job = Gaia.launch_job_async(f"""
                SELECT ra, dec, phot_g_mean_mag
                FROM gaiadr3.gaia_source
                WHERE ra  BETWEEN {ra_min}  AND {ra_max}
                AND   dec BETWEEN {dec_min} AND {dec_max}
                AND   phot_g_mean_mag IS NOT NULL
                AND   (classprob_dsc_combmod_star > 0.98 OR classprob_dsc_combmod_star IS NULL OR phot_g_mean_mag < 14)
            """)
            job.get_results().write(gaia_out, format="csv", overwrite=True)
            gaia = pd.read_csv(gaia_out).astype({"ra": float, "dec": float})
            print(f"  Gaia (fresh): {len(gaia)}")

        # kd-tree nearest-neighbor matching on the unit sphere
        faint_mask = gaia["phot_g_mean_mag"].values >= 11
        se_faint = catalog_psf[catalog_psf["MAG_CAL"] >= o2_hi].copy().reset_index(drop=True)

        dist, idx = cKDTree(_xyz(gaia["ra"].values[faint_mask], gaia["dec"].values[faint_mask])).query(
            _xyz(se_faint["RA"].values, se_faint["Dec"].values), k=1, workers=-1
        )
        sep = np.rad2deg(2 * np.arcsin(np.clip(dist / 2, -1, 1))) * 3600

        # Per-tile match-radius threshold from the separation histogram,
        # rather than a fixed arcsec cut
        counts, edges = np.histogram(sep[sep < 5], bins=200)
        centres = (edges[:-1] + edges[1:]) / 2
        pk = np.argmax(counts)
        below = counts[pk:] < max(3, 0.05 * counts[pk])
        thresh = float(centres[pk:][np.argmax(below)]) if below.any() else 2.0
        print(f"  Match threshold: {thresh:.3f}\"")

        within = sep < thresh
        matched_a = se_faint[within].copy().reset_index(drop=True)
        gaia_sub = gaia[faint_mask].reset_index(drop=True).iloc[idx[within]].reset_index(drop=True)
        matched_a["GAIA_RA"] = gaia_sub["ra"].values
        matched_a["GAIA_Dec"] = gaia_sub["dec"].values
        matched_a["GAIA_Gmag"] = gaia_sub["phot_g_mean_mag"].values
        matched_a["Sep_arcsec"] = sep[within]
        matched_a["TILE"] = tile_name
        matched_a["SOURCE"] = "Gaia_matched"
        print(f"  Gaia-matched: {len(matched_a)}")

        # Combine Gaia-matched + SE-only, dedupe, then check centering
        matched = (
            pd.concat([matched_a, se_outer2], ignore_index=True)
            .sort_values("Sep_arcsec", na_position="last")
            .drop_duplicates("NUMBER", keep="first")
            .reset_index(drop=True)
        )

        saturated = (matched["FLAGS"] & 4).astype(bool)
        image_data = hdul[0].data
        matched["CENTRED"] = [
            True if saturated.iloc[i] else _is_centered(image_data, r.X_IMAGE, r.Y_IMAGE)
            for i, r in enumerate(matched.itertuples())
        ]
        del image_data
        matched = matched[matched["CENTRED"]].drop(columns="CENTRED").reset_index(drop=True)

    matched["REGION"] = matched.apply(lambda r: get_region(r["GAIA_Gmag"], r["MAG_CAL"], mag_cuts), axis=1)

    for region, (lo, hi) in mag_cuts.items():
        sub = matched[matched["REGION"] == region]
        print(f"  {region:<15} {lo}\u2013{hi}: {len(sub):4d}" + (" (SE-only)" if region == "Outer_2" else ""))

    return fits_path, matched


def build_master_catalog(tile_names=None, mag_cuts=None, zp=None, use_cache=True, out_path=None):
    """
    Cross-match every tile, save per-tile CSVs, write the combined
    master catalog, and print a per-region summary.

    Returns (combined_df, all_matched) where all_matched is a list of
    (tile_name, fits_path, matched_df) tuples — the fits_path/matched_df
    pairs are what plot_source_overlay() needs.
    """
    tile_names = tile_names or config.TILE_NAMES
    mag_cuts = mag_cuts or config.MAG_CUTS
    out_path = out_path or config.MATCHED_CATALOG_PATH

    all_matched = []
    for tile_name in tile_names:
        fits_path, matched = cross_match_tile(tile_name, mag_cuts=mag_cuts, zp=zp, use_cache=use_cache)
        all_matched.append((tile_name, fits_path, matched))

    for tile_name, _, matched in all_matched:
        tile_csv = os.path.join(config.GAIA_OUTPUT_DIR, f"{tile_name}.csv")
        os.makedirs(config.GAIA_OUTPUT_DIR, exist_ok=True)
        matched.to_csv(tile_csv, index=False, float_format="%.4f")
        print(f"Saved: {tile_csv}  ({len(matched)} stars)")

    combined = pd.concat([m for _, _, m in all_matched], ignore_index=True)
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    combined.to_csv(out_path, index=False, float_format="%.4f")

    print(f"\nTotal: {len(combined)} stars across {len(tile_names)} tiles")
    for region, (lo, hi) in mag_cuts.items():
        sub = combined[combined["REGION"] == region]
        print(f"  {region:<15} {lo}\u2013{hi}: {len(sub):5d}" + (" (SE-only)" if region == "Outer_2" else ""))

    return combined, all_matched


def plot_source_overlay(all_matched, mag_cuts=None, colors=None, sizes=None, step=4, n_cols=3):
    """
    Grid of tile images (one per tile) with matched sources overlaid,
    color- and size-coded by region.

    `all_matched` is the list of (tile_name, fits_path, matched_df)
    tuples returned by build_master_catalog().
    """
    mag_cuts = mag_cuts or config.MAG_CUTS
    colors = colors or config.REGION_COLORS
    sizes = sizes or config.REGION_MARKER_SIZES

    n_rows = math.ceil(len(all_matched) / n_cols)
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(6 * n_cols, 6 * n_rows), facecolor="black")
    axes = np.array(axes).flatten()

    i = 0
    for i, (tile_name, fits_path, matched) in enumerate(all_matched):
        with fits.open(fits_path, memmap=True) as hdul:
            data = hdul[0].data[::step, ::step].copy()

        ax = axes[i]
        ax.set_facecolor("black")
        ax.imshow(data, origin="lower", cmap="Greys", vmin=np.percentile(data, 1), vmax=np.percentile(data, 99))
        ax.set_title(tile_name, color="white", fontsize=9)
        ax.axis("off")
        del data

        for region, (lo, hi) in mag_cuts.items():
            sub = matched[matched["REGION"] == region]
            if not len(sub):
                continue
            label = f"{region} G{lo}\u2013{hi} n={len(sub)}" + (" (SE)" if region == "Outer_2" else "")
            ax.scatter(
                sub["X_IMAGE"].values / step, sub["Y_IMAGE"].values / step,
                s=sizes[region], facecolors="none", edgecolors=colors[region],
                linewidths=0.9, alpha=0.85, label=label,
            )

        ax.legend(loc="upper right", fontsize=7, facecolor="#111", edgecolor="white", labelcolor="white")

    for j in range(i + 1, len(axes)):
        axes[j].set_visible(False)

    plt.suptitle(f"KiDS tiles \u2014 source overlay ({len(all_matched)} tiles)", color="white", fontsize=13)
    plt.tight_layout()
    finalize_plot(fig, "source_overlay.png", subdir="gaia_validation")


def main():
    combined, all_matched = build_master_catalog()
    plot_source_overlay(all_matched)


if __name__ == "__main__":
    main()