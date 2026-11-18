"""
Step: Point Source Identification.

Overlays SExtractor-detected point sources (CLASS_STAR above a fixed or
per-tile adaptive threshold) on each tile, and plots FLUX_RADIUS vs
MAG_AUTO to visualize the stellar locus.
"""

import os

import numpy as np
import pandas as pd
from astropy.io import fits
import matplotlib.pyplot as plt

from . import config
from .helpers import finalize_plot, setup_plot_mode


def _read_sextractor_catalog(cat_file):
    """Parse a SExtractor ASCII_HEAD .cat file into a DataFrame using its # header lines."""
    columns = []
    with open(cat_file, "r") as f:
        for line in f:
            if line.startswith("#"):
                parts = line.replace("#", "").strip().split()
                if len(parts) >= 2:
                    columns.append(parts[1])
            else:
                break
    return pd.read_csv(cat_file, comment="#", sep=r"\s+", names=columns)


def plot_point_sources(tile_names=None, class_star_thresh=0.84, downsample=4):
    """Overlay CLASS_STAR-selected point sources on each tile's downsampled image."""
    tile_names = tile_names or config.TILE_NAMES

    fig, axes = plt.subplots(1, len(tile_names), figsize=(6 * len(tile_names), 8),
                              facecolor="black", squeeze=False)
    axes = axes.flatten()

    for i, name in enumerate(tile_names):
        fpath = os.path.join(config.FITS_DIR, f"{name}.fits")
        cat_file = os.path.join(config.SE_IO_DIR, f"outparams_{name}.cat")

        catalog = _read_sextractor_catalog(cat_file)
        point_sources = catalog[catalog["CLASS_STAR"] >= class_star_thresh]
        print(f"{name}: {len(point_sources)} point sources")

        with fits.open(fpath, memmap=True) as hdul:
            display = hdul[0].data[::downsample, ::downsample].copy()

        xs = point_sources["X_IMAGE"].values / downsample
        ys = point_sources["Y_IMAGE"].values / downsample

        ax = axes[i]
        ax.set_facecolor("black")
        ax.imshow(display, origin="lower", cmap="Greys", vmin=np.percentile(display, 1), vmax=np.percentile(display, 99))
        ax.scatter(xs, ys, s=10, edgecolor="red", facecolors="none", linewidth=0.6, alpha=0.6)
        ax.set_title(f"{name}\n{len(point_sources)} point sources", color="white", fontsize=9)
        ax.axis("off")

    plt.suptitle(f"Point sources (CLASS_STAR \u2265 {class_star_thresh}) across KiDS tiles", color="white", fontsize=12)
    plt.tight_layout()
    finalize_plot(fig, "point_sources.png", subdir="point_sources")


def adaptive_class_star_threshold(counts, centres, lo=0.3, hi=0.9, clip=(0.5, 0.9)):
    """
    Find the CLASS_STAR value at the histogram valley between `lo` and `hi`
    (the natural galaxy/star separation point), clipped to `clip`.
    """
    valley = (centres > lo) & (centres < hi)
    return float(np.clip(centres[valley][np.argmin(counts[valley])], *clip))


def flux_radius_vs_mag(tile_names=None):
    """
    Per-tile adaptive CLASS_STAR threshold, concatenated catalog/point-source
    tables across all tiles, and a FLUX_RADIUS vs MAG_AUTO scatter plot.

    Returns (catalog_all, stars_all) DataFrames for use by later steps.
    """
    tile_names = tile_names or config.TILE_NAMES

    all_cats, all_stars = [], []

    for tile_name in tile_names:
        cat_file = os.path.join(config.SE_IO_DIR, f"outparams_{tile_name}.cat")
        catalog = _read_sextractor_catalog(cat_file)
        catalog["TILE"] = tile_name

        cs_counts, cs_edges = np.histogram(catalog["CLASS_STAR"], bins=100)
        cs_centres = (cs_edges[:-1] + cs_edges[1:]) / 2
        cs_thresh = adaptive_class_star_threshold(cs_counts, cs_centres)

        point_sources = catalog[catalog["CLASS_STAR"] >= cs_thresh]
        all_cats.append(catalog)
        all_stars.append(point_sources)
        print(f"{tile_name}: {len(catalog)} sources, {len(point_sources)} point sources (CS\u2265{cs_thresh:.3f})")

    catalog_all = pd.concat(all_cats, ignore_index=True)
    stars_all = pd.concat(all_stars, ignore_index=True)

    fig, ax = plt.subplots(figsize=(10, 7), facecolor="black")
    ax.set_facecolor("#0a0a0a")

    ax.scatter(catalog_all["FLUX_RADIUS"], catalog_all["MAG_AUTO"], s=8, color="white", alpha=0.15, label="All sources", zorder=1)
    ax.scatter(stars_all["FLUX_RADIUS"], stars_all["MAG_AUTO"], s=20, edgecolors="teal", facecolors="none",
               alpha=0.7, linewidths=0.6, label="Point sources", zorder=2)

    ax.set_xlim(-3, 40)
    ax.set_ylim(27, 8)
    ax.set_xlabel("FLUX_RADIUS (pixels)", color="white", fontsize=11)
    ax.set_ylabel("MAG_AUTO", color="white", fontsize=11)
    ax.set_title("FLUX_RADIUS vs MAG_AUTO", color="white", fontsize=12, pad=10)

    ax.tick_params(colors="white")
    for sp in ax.spines.values():
        sp.set_edgecolor("#444")
    ax.grid(True, color="white", alpha=0.07, linestyle="--")
    ax.set_axisbelow(True)
    ax.legend(fontsize=9, facecolor="#111", edgecolor="white", labelcolor="white")

    plt.tight_layout()
    finalize_plot(fig, "flux_radius_vs_mag.png", subdir="point_sources")

    return catalog_all, stars_all


def main():
    plot_point_sources()
    flux_radius_vs_mag()


if __name__ == "__main__":
    main()