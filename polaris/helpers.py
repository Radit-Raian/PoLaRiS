"""
Helper functions shared across two or more PoLaRiS pipeline steps.

Kept separate from any single step's module so that, e.g., stacking.py
and cutouts.py both import the same `is_valid_cutout` implementation
instead of maintaining independent copies that can drift apart.
"""

import os
import sys

import numpy as np
from astropy.io import fits

from . import config

# Cached plot-mode choice, set explicitly via setup_plot_mode() once at
# the start of a driver script. Every step's finalize_plot() call in a
# single run then shares that one choice.
_SHOW_PLOTS = None


def setup_plot_mode(show=None):
    """
    Set, once, whether pipeline plots are shown interactively or saved
    to disk. Call this explicitly at the start of your driver script,
    before running any pipeline step:

        from polaris.helpers import setup_plot_mode
        setup_plot_mode()            # prompts interactively
        setup_plot_mode(show=False)  # or set it directly, no prompt

    If `show` is left as None, the user is prompted — unless stdin
    isn't a real terminal (e.g. a scheduled/non-interactive run), in
    which case it defaults to False rather than blocking on input().
    """
    global _SHOW_PLOTS
    if show is not None:
        _SHOW_PLOTS = bool(show)
    elif not sys.stdin.isatty():
        print("Non-interactive session detected — plots will be saved "
              f"under {config.PLOTS_DIR} instead of displayed.")
        _SHOW_PLOTS = False
    else:
        answer = input("Show intermediate plots interactively? [y/N]: ").strip().lower()
        _SHOW_PLOTS = answer in ("y", "yes")

    if not _SHOW_PLOTS:
        print(f"Plots will be saved under {config.PLOTS_DIR}.")
    return _SHOW_PLOTS


def reset_plot_mode():
    """Clear the cached choice from setup_plot_mode(), e.g. between test runs."""
    global _SHOW_PLOTS
    _SHOW_PLOTS = None


def _current_plot_mode():
    """
    Internal: return the cached setup_plot_mode() choice.

    Only run_pipeline.py calls setup_plot_mode() (once, up front) to ask
    the interactive y/N question. An individual step run standalone
    (`python -m polaris.<step>`) never triggers that prompt — it just
    silently saves its plots under config.PLOTS_DIR, since there's no
    "whole run" to ask about.
    """
    global _SHOW_PLOTS
    if _SHOW_PLOTS is None:
        _SHOW_PLOTS = False
    return _SHOW_PLOTS


def finalize_plot(fig, filename, subdir=None):
    """
    Show `fig` interactively, or save it under config.PLOTS_DIR (in an
    optional `subdir`) depending on the choice made via
    setup_plot_mode() — then close it either way so figures don't pile
    up across many calls (e.g. per-star previews).
    """
    import matplotlib.pyplot as plt

    if _current_plot_mode():
        plt.show()
    else:
        out_dir = os.path.join(config.PLOTS_DIR, subdir) if subdir else config.PLOTS_DIR
        os.makedirs(out_dir, exist_ok=True)
        out_path = os.path.join(out_dir, filename)
        fig.savefig(out_path, dpi=150, bbox_inches="tight", facecolor=fig.get_facecolor())
        print(f"  Saved plot: {out_path}")
    plt.close(fig)


def get_patch(flux_radius, mag, patches):
    """
    Return (cutout_size, folder_name) for the first patch definition
    in `patches` whose (flux_radius, mag) window contains this star.

    `patches` is a list of (fr_min, fr_max, mag_min, mag_max, size, folder)
    tuples, as defined in config.PATCHES.

    Returns (None, None) if the star doesn't fall in any patch.
    """
    for fr_min, fr_max, m_min, m_max, size, folder in patches:
        if fr_min <= flux_radius <= fr_max and m_min <= mag <= m_max:
            return size, folder
    return None, None


def is_valid_cutout(data, threshold):
    """
    Return False if more than `threshold` fraction of pixels in `data`
    are NaN or exactly zero (indicating an edge/empty cutout).
    """
    bad = np.sum(~np.isfinite(data) | (data == 0))
    return (bad / data.size) <= threshold


def make_radius_map(ny, nx):
    """Build a 2D array giving each pixel's distance (px) from the array center."""
    cy, cx = (ny - 1) / 2.0, (nx - 1) / 2.0
    yy, xx = np.ogrid[:ny, :nx]
    return np.sqrt((yy - cy) ** 2 + (xx - cx) ** 2)


def annulus_median(data, r_map, r_inner, r_outer):
    """Median of the finite, positive pixel values within an annulus."""
    mask = (r_map >= r_inner) & (r_map <= r_outer)
    vals = data[mask]
    vals = vals[np.isfinite(vals) & (vals > 0)]
    return float(np.nanmedian(vals)) if len(vals) > 0 else np.nan


def peak_offset(data):
    """Distance (px) of the brightest finite pixel from the array's geometric center."""
    finite = np.where(np.isfinite(data), data, -np.inf)
    peak_y, peak_x = np.unravel_index(np.argmax(finite), data.shape)
    cy, cx = data.shape[0] // 2, data.shape[1] // 2
    return np.sqrt((peak_x - cx) ** 2 + (peak_y - cy) ** 2)


def load_masked_cutout(cutout_path, mult_path):
    """
    Load a raw cutout and apply its corresponding multiplied-mask FITS
    file (if present), converting masked-out or exactly-zero pixels
    to NaN so downstream statistics ignore them.
    """
    cutout = fits.getdata(cutout_path).astype(np.float32)
    data = cutout.copy()
    if os.path.exists(mult_path):
        mult = fits.getdata(mult_path).astype(np.float32)
        data[(mult == 0) & (cutout != 0)] = np.nan
    data[cutout == 0] = np.nan
    return data


def plot_stack(median_stack, mean_stack, stack_name, n_accepted, expected_size, r_inner, r_outer):
    """Log-stretched side-by-side display of a patch's median and mean stacks."""
    import matplotlib.pyplot as plt

    def log_stretch(img):
        out = img.copy().astype(np.float64)
        out[out <= 0] = np.nan
        return np.log10(out)

    med_log = log_stretch(median_stack)
    mean_log = log_stretch(mean_stack)

    finite_vals = med_log[np.isfinite(med_log)]
    vmax = np.nanpercentile(finite_vals, 99.9)
    vmin = vmax - 6.5  # 6.5 dex range matches paper dynamic range

    fig, axes = plt.subplots(1, 2, figsize=(12, 5.5))
    fig.patch.set_facecolor("black")
    fig.suptitle(
        f"{stack_name}   N={n_accepted} stars   "
        f"cutout={expected_size}px   norm annulus={r_inner}\u2013{r_outer}px",
        color="white", fontsize=12, y=1.01,
    )

    for ax, (img, title) in zip(axes, [(med_log, "Median Stack"), (mean_log, "Mean Stack")]):
        ax.set_facecolor("black")
        im = ax.imshow(img, origin="lower", cmap="viridis", vmin=vmin, vmax=vmax, interpolation="nearest")
        ax.set_title(title, color="white", fontsize=11)
        ax.axis("off")

        cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label("log\u2081\u2080 (Normalised Flux)", color="white", fontsize=9)
        cbar.ax.yaxis.set_tick_params(color="white")
        plt.setp(cbar.ax.yaxis.get_ticklabels(), color="white")

        cy, cx = img.shape[0] / 2, img.shape[1] / 2
        for r in [r_inner, r_outer]:
            circle = plt.Circle((cx, cy), r, color="red", fill=False, linewidth=0.8, linestyle="--")
            ax.add_patch(circle)

    plt.tight_layout()
    finalize_plot(fig, f"stack_{stack_name}.png", subdir="stacking")