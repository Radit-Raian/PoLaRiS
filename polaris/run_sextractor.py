"""
Step: SExtractor Setup and Run.

Writes a per-tile SExtractor configuration file, runs SExtractor to
detect sources, and previews the resulting segmentation maps.
"""

import glob
import math
import os
import subprocess

import numpy as np
import pandas as pd
from astropy.io import fits
import matplotlib.pyplot as plt

from . import config
from .helpers import finalize_plot, setup_plot_mode


def _resolve_pixscale(hdr, default=config.DEFAULT_PIXSCALE):
    """
    Read the WCS pixel scale from a FITS header, arcsec/pixel.

    Falls back through CD1_1 -> CDELT1 -> `default` rather than raising
    a KeyError, since not every survey/tile header uses the same WCS
    keyword convention.
    """
    if "CD1_1" in hdr:
        return abs(hdr["CD1_1"]) * 3600
    if "CDELT1" in hdr:
        return abs(hdr["CDELT1"]) * 3600
    return default


def setup_directories(fits_dir=None):
    """Ensure the SExtractor input/output directory exists."""
    fits_dir = fits_dir or config.FITS_DIR
    os.makedirs(config.SE_IO_DIR, exist_ok=True)


def run_sextractor_on_tiles(fits_dir=None):
    """
    Loop over every FITS tile, write its SExtractor config template
    (using per-tile SEEING/SATUR/PIXSCALE read from the header), and
    run SExtractor.
    """
    fits_dir = fits_dir or config.FITS_DIR
    setup_directories(fits_dir)

    fits_pathlist = sorted(glob.glob(os.path.join(fits_dir, "*.fits")))
    fits_list = [os.path.basename(f) for f in fits_pathlist]

    env = os.environ.copy()

    for i, (image_file, filename) in enumerate(zip(fits_pathlist, fits_list), start=1):
        base_name = os.path.splitext(filename)[0]
        cat_out_name = os.path.join(config.SE_IO_DIR, f"outparams_{base_name}.cat")
        template_path = os.path.join(config.SE_IO_DIR, f"{base_name}.SEtemplate")

        with fits.open(image_file, memmap=True) as hdul:
            hdr = hdul[0].header
            seeing = hdr.get("PSF_RAD", 0.7)  # arcsec, tile-specific
            satur = hdr.get("STATMAX", 3.5e-8) * 0.95  # 95% of brightest pixel
            pixscale = _resolve_pixscale(hdr)

        print(f"\n[{i}/{len(fits_list)}] {base_name}")
        print(f"  SEEING={seeing:.4f}\"  SATUR={satur:.3e}  PIXSCALE={pixscale:.4f}\"/px")

        with open(template_path, "w") as f:
            f.write(f"""CATALOG_NAME     {cat_out_name}

CATALOG_TYPE     ASCII_HEAD
PARAMETERS_NAME  {config.OUT_PARAM_PATH}

DETECT_TYPE      CCD
DETECT_MINAREA   5
DETECT_THRESH    2.0
ANALYSIS_THRESH  2.0

FILTER           Y
FILTER_NAME      {config.CONV_FILTER}

DEBLEND_NTHRESH  32
DEBLEND_MINCONT  0.001

CLEAN            Y
CLEAN_PARAM      1.0

WEIGHT_TYPE      NONE
WEIGHT_IMAGE     -

FLAG_IMAGE       flag.fits
FLAG_TYPE        OR

PHOT_APERTURES   10
PHOT_AUTOPARAMS  2.5,3.5
PHOT_PETROPARAMS 2.0,3.5
PHOT_AUTOAPERS   0.0,0.0

SATUR_LEVEL      {satur}
SATUR_KEY        SATURATE

MAG_ZEROPOINT    0.0
MAG_GAMMA        4.0
GAIN             1.0
GAIN_KEY         GAIN
PIXEL_SCALE      {pixscale}

SEEING_FWHM      {seeing}
STARNNW_NAME     {config.NNW_PATH}
PSF_NAME         {config.PSF_PATH}

CHECKIMAGE_TYPE  SEGMENTATION,BACKGROUND,BACKGROUND_RMS
CHECKIMAGE_NAME  {config.SE_IO_DIR}/segmap_{base_name}.fits,{config.SE_IO_DIR}/bkg_{base_name}.fits,{config.SE_IO_DIR}/bkg_rms_{base_name}.fits

BACK_TYPE        AUTO
BACK_SIZE        128
BACK_FILTERSIZE  3

MEMORY_OBJSTACK  50000
MEMORY_PIXSTACK  3000000
MEMORY_BUFSIZE   8000

VERBOSE_TYPE     NORMAL
HEADER_SUFFIX    .head
WRITE_XML        N
XML_NAME         sex.xml
XSL_URL          file:///usr/local/share/sextractor/sextractor.xsl
""")

        result = subprocess.run(
            [config.SOURCE_EXTRACTOR, image_file, "-c", template_path],
            capture_output=True, text=True, env=env,
        )

        if result.returncode == 0:
            cat = pd.read_csv(
                cat_out_name, sep=r"\s+", comment="#",
                names=["NUMBER", "MAG_AUTO", "MAGERR_AUTO", "X_IMAGE", "Y_IMAGE",
                       "FLAGS", "ELLIPTICITY", "CLASS_STAR", "FLUX_RADIUS"],
                index_col=False,
            )
            n_stars = (cat["CLASS_STAR"] >= 0.84).sum()
            print(f"  \u2713 Done \u2014 {len(cat)} detections, {n_stars} stars (CLASS_STAR\u22650.84)")
        else:
            print(f"  \u2717 Failed (code {result.returncode})")
            print(result.stderr[-500:])

    print(f"\nSExtractor ran for all {len(fits_list)} FITS files.")


def preview_segmaps(fits_dir=None, downsample=4):
    """Grid preview of every tile's SExtractor segmentation map."""
    fits_dir = fits_dir or config.FITS_DIR
    segmap_files = sorted(glob.glob(os.path.join(config.SE_IO_DIR, "segmap_KIDS_*.fits")))

    print(f"Found {len(segmap_files)} segmaps:")
    for f in segmap_files:
        print(f"  {os.path.basename(f)}")

    if not segmap_files:
        print("No segmaps found — run run_sextractor_on_tiles() first.")
        return

    n_cols = min(3, len(segmap_files))
    n_rows = math.ceil(len(segmap_files) / n_cols)

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(6 * n_cols, 6 * n_rows), facecolor="black")
    axes = np.array(axes).flatten()

    i = 0
    for i, spath in enumerate(segmap_files):
        try:
            with fits.open(spath, memmap=True) as hdul:
                seg = hdul[0].data[::downsample, ::downsample].copy()
                name = os.path.basename(spath).replace("segmap_", "").replace(".fits", "")

            ax = axes[i]
            ax.set_facecolor("black")
            im = ax.imshow(seg, origin="lower", cmap="tab20")
            plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="Object ID")
            ax.set_title(name, color="white", fontsize=9)
            ax.axis("off")
        except Exception as e:
            axes[i].set_visible(False)
            print(f"Skipping {os.path.basename(spath)}: {e}")

    for j in range(i + 1, len(axes)):
        axes[j].set_visible(False)

    plt.suptitle(f"Segmentation maps ({len(segmap_files)} tiles)", color="white", fontsize=13)
    plt.tight_layout()
    finalize_plot(fig, "segmaps_preview.png", subdir="run_sextractor")


def main():
    run_sextractor_on_tiles()
    preview_segmaps()


if __name__ == "__main__":
    main()