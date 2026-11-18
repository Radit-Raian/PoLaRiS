"""
Step: Cutout Generation.

Extracts a fixed-size cutout around each selected star, sized according
to which patch (Core / Intermediate / Outer_1 / Outer_2) it falls into.

`extract_cutout()` is the single source of truth for turning a
(tile image, x, y, size) into a saved array — masking.py imports and
reuses it directly instead of re-implementing Cutout2D extraction, so
the two steps can't drift out of sync with each other.

Running this module standalone (`python -m polaris.cutouts`) is
intended only as a quick QA pass on cutout extraction before running
the heavier masking step, which regenerates cutouts alongside their
masks — see masking.py.
"""

import os
import random

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.nddata import Cutout2D
import matplotlib.pyplot as plt

from . import config
from .helpers import get_patch, is_valid_cutout, finalize_plot, setup_plot_mode


def extract_cutout(img, x, y, size, dtype="float32", fill_value=0):
    """
    Extract a (size x size) cutout centered at pixel (x, y) from `img`.

    Uses mode='partial' so stars near a tile edge still produce a
    correctly-shaped array, with `fill_value` in the region outside
    the tile.
    """
    cut = Cutout2D(img, (x, y), (size, size), mode="partial", fill_value=fill_value)
    return cut.data.astype(dtype)


def generate_cutouts(matched_catalog_path=None, tile_names=None, patches=None,
                      edge_margin_frac=0.55):
    """
    For every star in the matched catalog, determine its patch, extract
    its cutout, validate it, and save it under
    output/patches/<folder>/cutouts/ (see config.patch_dir).

    Returns a dict of per-tile (saved, skipped) counts.
    """
    matched_catalog_path = matched_catalog_path or config.MATCHED_CATALOG_PATH
    tile_names = tile_names or config.TILE_NAMES
    patches = patches or config.PATCHES

    for folder in config.UNIQUE_FOLDERS:
        os.makedirs(os.path.join(config.patch_dir(folder), "cutouts"), exist_ok=True)

    matched_catalog = pd.read_csv(matched_catalog_path)
    counts = {}

    for tile_name in tile_names:
        fits_path = os.path.join(config.FITS_DIR, f"{tile_name}.fits")
        with fits.open(fits_path, memmap=True) as hdul:
            img = hdul[0].data

        ny, nx = img.shape
        tile_stars = matched_catalog[matched_catalog["TILE"] == tile_name]

        saved = skipped = 0
        for _, row in tile_stars.iterrows():
            x, y = row["X_IMAGE"] - 1, row["Y_IMAGE"] - 1
            size, folder = get_patch(row["FLUX_RADIUS"], row["GAIA_Gmag"], patches)
            if folder is None:
                continue

            margin = size * edge_margin_frac
            if not (margin < x < nx - margin and margin < y < ny - margin):
                skipped += 1
                continue

            data = extract_cutout(img, x, y, size)
            if not is_valid_cutout(data, config.MAX_EMPTY_FRACTION):
                skipped += 1
                continue

            out_name = f"{tile_name}_{int(row['NUMBER'])}_{size}.fits"
            out_path = os.path.join(config.patch_dir(folder), "cutouts", out_name)
            fits.writeto(out_path, data, overwrite=True)
            saved += 1

        counts[tile_name] = (saved, skipped)
        print(f"{tile_name}: saved={saved}  skipped={skipped}")

    total_saved = sum(s for s, _ in counts.values())
    total_skipped = sum(k for _, k in counts.values())
    print(f"\nTotal: saved={total_saved}  skipped={total_skipped}")
    return counts


def preview_cutouts(n=5, seed=0):
    """Show `n` random saved cutouts per patch folder, percentile-scaled."""
    rng = random.Random(seed)

    for folder in config.UNIQUE_FOLDERS:
        cutout_dir = os.path.join(config.patch_dir(folder), "cutouts")
        if not os.path.isdir(cutout_dir):
            continue

        files = [f for f in os.listdir(cutout_dir) if f.endswith(".fits")]
        if not files:
            continue

        sample = rng.sample(files, min(n, len(files)))
        fig, axes = plt.subplots(1, len(sample), figsize=(4 * len(sample), 4))
        if len(sample) == 1:
            axes = [axes]

        for ax, fname in zip(axes, sample):
            data = fits.getdata(os.path.join(cutout_dir, fname))
            vmin, vmax = np.nanpercentile(data, [1, 99])
            ax.imshow(data, origin="lower", cmap="viridis", vmin=vmin, vmax=vmax)
            ax.set_title(fname, fontsize=7)
            ax.axis("off")

        fig.suptitle(folder)
        plt.tight_layout()
        finalize_plot(fig, f"cutout_preview_{folder}.png", subdir="cutouts")


def main():
    generate_cutouts()
    preview_cutouts()


if __name__ == "__main__":
    main()