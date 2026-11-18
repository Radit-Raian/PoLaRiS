"""
Step: Stacking.

Loads every masked cutout per patch, applies rejection criteria
(wrong shape, too few valid pixels, too heavily masked, miscentered,
bad normalization), applies a random continuous-angle rotation to
decorrelate non-circular PSF features (diffraction spikes, chip
artifacts) across stars, normalizes by an annulus median, and stacks.

Saves median stack, mean stack, and a per-pixel uncertainty map
(Garate et al.) as FITS files under output/patches/<patch>/stacked/.
"""

import os
import random

import numpy as np
import pandas as pd
from astropy.io import fits
from scipy.ndimage import rotate as ndi_rotate

from . import config
from .helpers import make_radius_map, annulus_median, peak_offset, load_masked_cutout, plot_stack, setup_plot_mode


def _rotate_random_angle(data):
    """
    Rotate `data` by a random angle in [0, 360). Continuous angles (as
    opposed to random 90-degree multiples) are needed because features
    aligned to the pixel grid \u2014 diffraction spikes, chip seams \u2014 are
    NOT decorrelated by 90-degree-only rotation.

    order=1 (bilinear) avoids ringing near the sharp edges of the NaN
    mask. NaNs propagate a ~1px halo into neighboring pixels through
    interpolation; cval=np.nan fills the corners rotated in from
    outside the original frame.
    """
    angle = random.uniform(0, 360)
    return ndi_rotate(data, angle, reshape=False, order=1, mode="constant", cval=np.nan)


def stack_patch(stack_name, patch_list, show_plots=True):
    """Run the full stacking pipeline for one patch (or group of patches)."""
    stacked_dir = os.path.join(config.patch_dir(stack_name), "stacked")
    os.makedirs(stacked_dir, exist_ok=True)

    expected_size = config.PATCH_SIZES[stack_name]
    r_inner, r_outer = config.NORM_ANNULUS[stack_name]
    max_masked = config.MAX_MASKED_FRAC[stack_name]
    centre_tol = config.CENTRE_TOL[stack_name]

    stack_list = []
    n_total = n_wrong_shape = n_low_valid = n_high_masked = 0
    n_bad_centre = n_bad_norm = n_accepted = 0

    for patch in patch_list:
        cutouts_dir = os.path.join(config.patch_dir(patch), "cutouts")
        multiplied_dir = os.path.join(config.patch_dir(patch), "multiplied")

        if not os.path.isdir(cutouts_dir):
            print(f"  [WARN] Missing: {cutouts_dir}")
            continue

        fnames = sorted(f for f in os.listdir(cutouts_dir) if f.endswith(".fits"))

        for fname in fnames:
            n_total += 1
            data = load_masked_cutout(
                os.path.join(cutouts_dir, fname),
                os.path.join(multiplied_dir, fname),
            )

            if data.shape != (expected_size, expected_size):
                n_wrong_shape += 1
                continue

            r_map = make_radius_map(*data.shape)
            valid_frac = np.sum(np.isfinite(data)) / data.size

            if valid_frac < config.MIN_VALID_FRACTION:
                n_low_valid += 1
                continue

            if (1.0 - valid_frac) > max_masked:
                n_high_masked += 1
                continue

            if peak_offset(data) > centre_tol:
                n_bad_centre += 1
                continue

            norm_val = annulus_median(data, r_map, r_inner, r_outer)
            if not (np.isfinite(norm_val) and norm_val > 0):
                n_bad_norm += 1
                continue

            data = _rotate_random_angle(data)

            stack_list.append(data / norm_val)
            n_accepted += 1

    summary = dict(
        Stack=stack_name, Total=n_total, Accepted=n_accepted,
        Rej_shape=n_wrong_shape, Rej_valid=n_low_valid,
        Rej_masked=n_high_masked, Rej_centre=n_bad_centre,
        Rej_norm=n_bad_norm,
    )

    if n_accepted == 0:
        print(f"[SKIP] {stack_name}: no valid PSFs \u2014 check summary below")
        return summary

    stack_array = np.array(stack_list)
    median_stack = np.nanmedian(stack_array, axis=0)
    mean_stack = np.nanmean(stack_array, axis=0)

    q16 = np.nanpercentile(stack_array, 16, axis=0)
    q84 = np.nanpercentile(stack_array, 84, axis=0)
    n_contrib = np.sum(np.isfinite(stack_array), axis=0).astype(float)
    n_contrib[n_contrib == 0] = np.nan
    sigma_stack = ((q84 - q16) / 2.0) * 1.253 / np.sqrt(n_contrib)

    for arr, tag in [(median_stack, "median"), (mean_stack, "mean"), (sigma_stack, "sigma")]:
        fits.writeto(
            os.path.join(stacked_dir, f"{tag}_psf_{stack_name}.fits"),
            arr.astype("float32"), overwrite=True,
        )

    if show_plots:
        plot_stack(median_stack, mean_stack, stack_name, n_accepted, expected_size, r_inner, r_outer)

    return summary


def run_stacking(stack_groups=None, show_plots=True):
    """Run stack_patch() for every group in `stack_groups` and print a summary table."""
    stack_groups = stack_groups or config.STACK_GROUPS
    summary_rows = [stack_patch(name, patches, show_plots=show_plots) for name, patches in stack_groups.items()]

    print("\n" + "=" * 70)
    print(pd.DataFrame(summary_rows).to_string(index=False))
    return summary_rows


def main():
    run_stacking()


if __name__ == "__main__":
    main()