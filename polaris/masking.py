"""
Step: Masking.

For every star, extracts the science cutout and the matching segmap
cutout (reusing cutouts.extract_cutout so this step never duplicates
the Cutout2D logic), identifies the target object at the cutout
center, builds a mask that zeroes out contaminating neighboring
objects, and saves all five per-star products: cutout, segmap cutout,
primary mask, final mask, and the multiplied (masked) image.

This step regenerates cutouts from scratch (wiping any output left by
cutouts.generate_cutouts) because the segmap cutout must be extracted
with the exact same center/size as the science cutout, and centering
can shift slightly once the target-object check below is applied.
Running cutouts.generate_cutouts() first is optional QA, not a
prerequisite.
"""

import os
import shutil

import numpy as np
import pandas as pd
from astropy.io import fits
from scipy.ndimage import binary_dilation, generate_binary_structure
import matplotlib.pyplot as plt

from . import config
from .cutouts import extract_cutout
from .helpers import get_patch, is_valid_cutout, finalize_plot, setup_plot_mode


def _wipe_and_recreate(patch_root):
    """Wipe and recreate each MASK_SUBDIRS folder under `patch_root` (a full path from config.patch_dir)."""
    for sub in config.MASK_SUBDIRS:
        path = os.path.join(patch_root, sub)
        if os.path.isdir(path):
            shutil.rmtree(path)
        os.makedirs(path, exist_ok=True)


def build_masks(matched_catalog_path=None, tile_names=None, patches=None,
                 edge_margin_frac=0.55, centre_tol_px=5):
    """
    For every star: extract cutout + segmap cutout, build primary/final
    masks, compute the multiplied image, and save all five products.

    Returns a dict of per-patch accepted counts.
    """
    matched_catalog_path = matched_catalog_path or config.MATCHED_CATALOG_PATH
    tile_names = tile_names or config.TILE_NAMES
    patches = patches or config.PATCHES

    for folder in config.UNIQUE_FOLDERS:
        _wipe_and_recreate(config.patch_dir(folder))

    matched_catalog = pd.read_csv(matched_catalog_path)
    struct = generate_binary_structure(2, 2)  # full 8-connectivity
    patch_counts = {folder: 0 for folder in config.UNIQUE_FOLDERS}

    for tile_name in tile_names:
        fits_path = os.path.join(config.FITS_DIR, f"{tile_name}.fits")
        segmap_path = os.path.join(config.SE_IO_DIR, f"segmap_{tile_name}.fits")

        with fits.open(fits_path, memmap=True) as hdul:
            img = hdul[0].data
        with fits.open(segmap_path, memmap=True) as hdul:
            seg = hdul[0].data

        ny, nx = img.shape
        tile_stars = matched_catalog[matched_catalog["TILE"] == tile_name]

        for _, row in tile_stars.iterrows():
            x, y = row["X_IMAGE"] - 1, row["Y_IMAGE"] - 1
            size, folder = get_patch(row["FLUX_RADIUS"], row["GAIA_Gmag"], patches)
            if folder is None:
                continue

            margin = size * edge_margin_frac
            if not (margin < x < nx - margin and margin < y < ny - margin):
                continue

            cutout_data = extract_cutout(img, x, y, size, dtype="float32")
            if not is_valid_cutout(cutout_data, config.MAX_EMPTY_FRACTION):
                continue

            seg_data = extract_cutout(seg, x, y, size, dtype="int32", fill_value=0)

            cy, cx = size // 2, size // 2
            target_id = seg_data[cy, cx]
            if target_id == 0:
                continue  # center pixel isn't on any detected object — skip

            star_id = f"{tile_name}_{int(row['NUMBER'])}_{size}"

            # Primary mask: star + background = 0, everything else = 1
            primary_mask = np.ones_like(seg_data, dtype=np.uint8)
            primary_mask[(seg_data == 0) | (seg_data == target_id)] = 0

            # Grow contamination flags outward so faint wings of nearby
            # objects are also excluded, not just their detected core
            contamination = (seg_data != 0) & (seg_data != target_id)
            contamination = binary_dilation(contamination, structure=struct, iterations=config.DILATION_ITERATIONS)
            primary_mask[contamination] = 1

            final_mask = 1 - primary_mask  # inverse: 1 = keep, 0 = mask out
            multiplied = cutout_data * final_mask

            for arr, subdir in [
                (cutout_data, "cutouts"),
                (seg_data, "segmaps"),
                (primary_mask, "primary_masking"),
                (final_mask, "final_masking"),
                (multiplied, "multiplied"),
            ]:
                out_path = os.path.join(config.patch_dir(folder), subdir, f"{star_id}.fits")
                fits.writeto(out_path, arr, overwrite=True)

            patch_counts[folder] += 1

    for folder, count in patch_counts.items():
        print(f"{folder}: {count} accepted")

    return patch_counts


def preview_masking(n=5, seed=0):
    """5-panel preview (cutout, segmap, primary mask, final mask, multiplied) for random stars per patch."""
    import random
    rng = random.Random(seed)

    for folder in config.UNIQUE_FOLDERS:
        cutout_dir = os.path.join(config.patch_dir(folder), "cutouts")
        if not os.path.isdir(cutout_dir):
            continue

        files = [f for f in os.listdir(cutout_dir) if f.endswith(".fits")]
        if not files:
            continue

        sample = rng.sample(files, min(n, len(files)))

        for fname in sample:
            fig, axes = plt.subplots(1, 5, figsize=(20, 4))
            panels = [
                ("cutouts", "Cutout", "viridis"),
                ("segmaps", "Segmap", "tab20"),
                ("primary_masking", "Primary Mask", "gray"),
                ("final_masking", "Final Mask", "gray"),
                ("multiplied", "Masked Image", "viridis"),
            ]
            for ax, (subdir, title, cmap) in zip(axes, panels):
                data = fits.getdata(os.path.join(config.patch_dir(folder), subdir, fname))
                if cmap == "viridis":
                    vmin, vmax = np.nanpercentile(data[data != 0], [1, 99]) if np.any(data != 0) else (0, 1)
                    ax.imshow(data, origin="lower", cmap=cmap, vmin=vmin, vmax=vmax)
                else:
                    ax.imshow(data, origin="lower", cmap=cmap)
                ax.set_title(title, fontsize=9)
                ax.axis("off")

            fig.suptitle(f"{folder} \u2014 {fname}")
            plt.tight_layout()
            star_id = os.path.splitext(fname)[0]
            finalize_plot(fig, f"mask_preview_{folder}_{star_id}.png", subdir="masking")


def main():
    build_masks()
    preview_masking()


if __name__ == "__main__":
    main()