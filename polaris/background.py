"""
Step: Background Subtraction.

Visualizes raw / background / background-RMS / background-subtracted
images for each tile, using the SExtractor CHECKIMAGE outputs produced
by run_sextractor.py.
"""

import os

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

from . import config
from .helpers import finalize_plot, setup_plot_mode


def plot_background_subtraction(tile_names=None, downsample=4):
    """Grid of raw / background / RMS / subtracted panels, one row per tile."""
    tile_names = tile_names or config.TILE_NAMES

    fig, axes = plt.subplots(len(tile_names), 4, figsize=(24, 6 * len(tile_names)),
                              facecolor="black", squeeze=False)

    col_titles = ["Original", "Background", "Background RMS", "Background Subtracted"]

    for i, name in enumerate(tile_names):
        fpath = os.path.join(config.FITS_DIR, f"{name}.fits")
        bkg_path = os.path.join(config.SE_IO_DIR, f"bkg_{name}.fits")
        rms_path = os.path.join(config.SE_IO_DIR, f"bkg_rms_{name}.fits")

        with fits.open(fpath, memmap=True) as h1, \
             fits.open(bkg_path, memmap=True) as h2, \
             fits.open(rms_path, memmap=True) as h3:
            raw = h1[0].data[::downsample, ::downsample].copy()
            bkg = h2[0].data[::downsample, ::downsample].copy()
            rms = h3[0].data[::downsample, ::downsample].copy()

        sub = raw - bkg

        for j, (data, title) in enumerate(zip([raw, bkg, rms, sub], col_titles)):
            ax = axes[i, j]
            ax.set_facecolor("black")
            ax.imshow(data, origin="lower", cmap="Greys", vmin=np.percentile(data, 1), vmax=np.percentile(data, 99))
            if i == 0:
                ax.set_title(title, color="white", fontsize=11)
            if j == 0:
                ax.set_ylabel(name, color="white", fontsize=9)
            ax.axis("off")

    plt.suptitle("KiDS tiles \u2014 raw / background / RMS / subtracted", color="white", fontsize=13)
    plt.tight_layout()
    finalize_plot(fig, "background_subtraction.png", subdir="background")


def main():
    plot_background_subtraction()


if __name__ == "__main__":
    main()