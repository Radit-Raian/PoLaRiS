"""
Step: Data Preprocessing.

Visualizes the raw FITS tiles (downsampled for quick inspection) so data
quality issues are caught before any downstream processing runs.
"""

import glob
import math
import os

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

from . import config
from .helpers import finalize_plot, setup_plot_mode


def inspect_raw_fits(fits_dir=None, downsample=4):
    """Plot a downsampled grayscale preview of every FITS tile in `fits_dir`."""
    fits_dir = fits_dir or config.FITS_DIR
    fits_files = sorted(glob.glob(os.path.join(fits_dir, "KIDS_*.fits")))

    print(f"Found {len(fits_files)} tiles:")
    for f in fits_files:
        print(f"  {os.path.basename(f)}")

    if not fits_files:
        print("No tiles found — check config.FITS_DIR.")
        return

    n_cols = 4
    n_rows = math.ceil(len(fits_files) / n_cols)

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(6 * n_cols, 6 * n_rows), facecolor="black")
    axes = np.array(axes).flatten()

    i = 0
    for i, fpath in enumerate(fits_files):
        with fits.open(fpath, memmap=True) as hdul:
            data = hdul[0].data[::downsample, ::downsample].copy()
            name = hdul[0].header.get("OBJECT", os.path.basename(fpath))

        ax = axes[i]
        ax.set_facecolor("black")
        ax.imshow(data, origin="lower", cmap="Greys", vmin=np.percentile(data, 1), vmax=np.percentile(data, 99))
        ax.set_title(name, color="white", fontsize=9)
        ax.axis("off")

    for j in range(i + 1, len(axes)):
        axes[j].set_visible(False)

    plt.suptitle(f"KiDS tiles ({len(fits_files)} total)", color="white", fontsize=13)
    plt.tight_layout()
    finalize_plot(fig, "raw_fits_preview.png", subdir="preprocessing")


def main():
    inspect_raw_fits()


if __name__ == "__main__":
    main()