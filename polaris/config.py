"""
Central configuration for the PoLaRiS pipeline.

Every path, per-patch parameter, and plotting constant used by the other
modules is defined here once, so a single edit (e.g. changing a survey's
data directory, or widening a patch's magnitude cut) is guaranteed to be
picked up consistently by every pipeline step.

Directory layout is auto-detected relative to this file, so cloning the
repo anywhere and dropping data into input/ works with zero path edits:

    <project root>/
    ├── input/
    │   ├── FITS_images/   survey tiles, e.g. KIDS_9.3_-31.2.fits (FITS_DIR)
    │   └── sextractor/    gauss_*.conv, default.nnw/.psf/.param
    ├── output/
    │   ├── se_in_out/         SExtractor catalogs, segmaps, backgrounds
    │   ├── patches/<name>/    cutouts / masks / stacked, per patch
    │   ├── gaia_cache/<tile>/ Gaia DR3 query results, auto-fetched and
    │   │                      cached per tile by gaia_validation.py —
    │   │                      not something you provide yourself
    │   └── pipeline_plots/<step>/   saved figures (see helpers.finalize_plot)
    └── polaris/           this package

Run `python -m polaris.config` once to scaffold empty input/output
folders if they don't exist yet. Every path below can still be
overridden with an environment variable for one-off runs.
"""

import os
from pathlib import Path

# ─────────────────────────────────────────────────────────────────────────
# Project layout — auto-detected, no hardcoded paths.
# PROJECT_ROOT is the folder that contains this `polaris` package (one
# level up from this file).
# ─────────────────────────────────────────────────────────────────────────
PROJECT_ROOT = Path(__file__).resolve().parent.parent

INPUT_DIR = Path(os.environ.get("POLARIS_INPUT_DIR", PROJECT_ROOT / "input"))
OUTPUT_DIR = Path(os.environ.get("POLARIS_OUTPUT_DIR", PROJECT_ROOT / "output"))

# ─────────────────────────────────────────────────────────────────────────
# Input — put your data here; nothing else needs to change.
# ─────────────────────────────────────────────────────────────────────────
FITS_DIR = str(Path(os.environ.get("POLARIS_FITS_DIR", INPUT_DIR / "FITS_images")))
SEXTRACTOR_INPUT_DIR = Path(os.environ.get("POLARIS_SEXTRACTOR_INPUT_DIR", INPUT_DIR / "sextractor"))

# SExtractor binary (still needs to be installed / on PATH) and the
# static support files you place in input/sextractor/.
SOURCE_EXTRACTOR = os.environ.get("POLARIS_SEX_BIN", "/usr/local/bin/sex")
CONV_FILTER = str(SEXTRACTOR_INPUT_DIR / "gauss_3.0_5x5.conv")
NNW_PATH = str(SEXTRACTOR_INPUT_DIR / "default.nnw")
PSF_PATH = str(SEXTRACTOR_INPUT_DIR / "default.psf")
OUT_PARAM_PATH = str(SEXTRACTOR_INPUT_DIR / "default.param")

# ─────────────────────────────────────────────────────────────────────────
# Output — every pipeline product lands under here, organized by step,
# so nothing is written loose into the project root anymore.
# ─────────────────────────────────────────────────────────────────────────
SE_IO_DIR = str(OUTPUT_DIR / "se_in_out")
PATCHES_DIR = str(OUTPUT_DIR / "patches")

# Where step outputs go when the user opts out of interactive plots
# (see helpers.setup_plot_mode / helpers.finalize_plot). Each step saves
# its figures into a named subfolder under here instead of popping up
# a window.
PLOTS_DIR = str(Path(os.environ.get("POLARIS_PLOTS_DIR", OUTPUT_DIR / "pipeline_plots")))

# Where gaia_validation.py writes its combined/per-tile matched catalogs
# and (via PLOTS_DIR) its overlay plot — kept together under one folder.
GAIA_OUTPUT_DIR = str(OUTPUT_DIR / "gaia_validation")

# Per-tile cache of raw Gaia DR3 query results. This is pipeline-generated
# cache, not something you supply — gaia_validation.py queries Gaia once
# per tile and writes the result here so repeat runs don't re-query. It
# lives under output/, not input/, for that reason.
GAIA_CACHE_DIR = str(OUTPUT_DIR / "gaia_cache")


def patch_dir(folder_name: str) -> str:
    """Root output directory for one patch's products (cutouts/segmaps/masks/stacked)."""
    return os.path.join(PATCHES_DIR, folder_name)


def scaffold_project():
    """
    Create the expected input/output folder skeleton (empty, with
    .gitkeep placeholders) if it doesn't exist yet. Safe to re-run.
    Run via: python -m polaris.config
    """
    input_dirs = [
        Path(FITS_DIR),
        SEXTRACTOR_INPUT_DIR,
    ]
    output_dirs = [Path(SE_IO_DIR), Path(GAIA_OUTPUT_DIR), Path(GAIA_CACHE_DIR)]
    for folder in UNIQUE_FOLDERS:
        output_dirs.append(Path(patch_dir(folder)))
        for sub in MASK_SUBDIRS + ["stacked"]:
            output_dirs.append(Path(patch_dir(folder)) / sub)
    for step in ["preprocessing", "run_sextractor", "background", "point_sources",
                 "star_selection", "gaia_validation", "cutouts", "masking",
                 "stacking", "radial_profile"]:
        output_dirs.append(Path(PLOTS_DIR) / step)

    for d in input_dirs + output_dirs:
        d.mkdir(parents=True, exist_ok=True)
        gitkeep = d / ".gitkeep"
        if not any(d.iterdir()):
            gitkeep.touch()

    print(f"Input directories ready under:  {INPUT_DIR}")
    for d in input_dirs:
        print(f"  {d}")
    print(f"\nOutput directories ready under: {OUTPUT_DIR}")
    print("  (patches/, se_in_out/, gaia_validation/, gaia_cache/, pipeline_plots/ — created on demand)")
    print(f"\nDrop your survey tiles into {FITS_DIR}")
    print(f"and your SExtractor support files (conv filter, .nnw, .psf, .param) into {SEXTRACTOR_INPUT_DIR}.")
    print("Gaia results are fetched and cached automatically — nothing to provide there.")


# ─────────────────────────────────────────────────────────────────────────
# Survey tiles
# ─────────────────────────────────────────────────────────────────────────
TILE_NAMES = [
    "KIDS_9.3_-31.2",
]

# Fallback pixel scale (arcsec/pixel) used only if a tile's FITS header
# lacks both CD1_1 and CDELT1 WCS keywords.
DEFAULT_PIXSCALE = 0.2

# Gaia cross-match catalogue produced by gaia_validation.py
MATCHED_CATALOG_PATH = os.path.join(GAIA_OUTPUT_DIR, "master.csv")


def tile_gaia_dir(tile_name: str) -> str:
    """Per-tile directory for cached Gaia query results (output cache, not input)."""
    return os.path.join(GAIA_CACHE_DIR, f"{tile_name}/")


# ─────────────────────────────────────────────────────────────────────────
# Gaia validation
#
# Magnitude-only region assignment (Core/Intermediate/Outer_1/Outer_2),
# used to label every matched star's REGION column. This is the
# authoritative patch assignment for a star — NOT the flux-radius-based
# PATCHES windows below, which star_selection.py / cutouts.py currently
# use instead. These two schemes can disagree on individual stars; see
# the note in star_selection.py.
# ─────────────────────────────────────────────────────────────────────────
GAIA_ZP = -0.044  # zero-point applied to SExtractor MAG_AUTO -> MAG_CAL

MAG_CUTS = {
    "Core": (16, 21),
    "Intermediate": (14, 16),
    "Outer_1": (11, 14),
    "Outer_2": (8, 11),
}

REGION_COLORS = {
    "Core": "teal",
    "Intermediate": "yellow",
    "Outer_1": "orange",
    "Outer_2": "red",
}

REGION_MARKER_SIZES = {
    "Core": 30,
    "Intermediate": 60,
    "Outer_1": 100,
    "Outer_2": 180,
}

# ─────────────────────────────────────────────────────────────────────────
# Star-selection patches
#
# Each tuple is (flux_radius_min, flux_radius_max, mag_min, mag_max,
# cutout_size_px, folder_name). A star is assigned to the first patch
# whose flux-radius/magnitude window contains it.
#
# ⚠ UNRECONCILED WITH gaia_validation.py: gaia_validation.py assigns
# each star's definitive REGION using MAG_CUTS (magnitude only, no
# flux radius) via get_region(). cutouts.py currently still calls
# get_patch() against PATCHES below, independently of that REGION
# column. Until these are reconciled, a star's folder here may not
# match its REGION in master.csv. See star_selection.py for the
# not-yet-resolved TODO.
# ─────────────────────────────────────────────────────────────────────────
PATCHES = [
    (0.0, 3.25, 16.0, 18.0, 200, "Core"),
    (0.0, 4.0, 14.5, 15.0, 500, "Intermediate"),
    (2.0, 8.0, 12.5, 14.5, 500, "Intermediate"),
    (8.0, 25.0, 11.0, 13.25, 1000, "Outer_1"),
    (25.0, 175.0, 8.0, 11.0, 2000, "Outer_2"),
]

UNIQUE_FOLDERS = list(dict.fromkeys(p[5] for p in PATCHES))

# Cutout size in pixels per patch (derived from PATCHES, kept explicit
# for readability in stacking.py / radial_profile.py)
PATCH_SIZES = {
    "Core": 200,
    "Intermediate": 500,
    "Outer_1": 1000,
    "Outer_2": 2000,
}

# Maximum fraction of a cutout that may be NaN/zero before it's rejected
MAX_EMPTY_FRACTION = 0.20

# ─────────────────────────────────────────────────────────────────────────
# Masking
# ─────────────────────────────────────────────────────────────────────────
MASK_SUBDIRS = ["cutouts", "segmaps", "primary_masking", "final_masking", "multiplied"]
DILATION_ITERATIONS = 4  # binary dilation of the contaminating-object mask

# ─────────────────────────────────────────────────────────────────────────
# Stacking
# ─────────────────────────────────────────────────────────────────────────
# Normalization annulus (inner, outer radius in px) per patch
NORM_ANNULUS = {
    "Core": (10, 15),
    "Intermediate": (100, 140),
    "Outer_1": (200, 220),
    "Outer_2": (400, 410),
}

# Maximum allowed masked-pixel fraction per patch, at stacking time
MAX_MASKED_FRAC = {
    "Core": 0.20,
    "Intermediate": 0.30,
    "Outer_1": 0.50,
    "Outer_2": 0.80,
}

# Maximum allowed offset (px) of a cutout's brightest pixel from its
# geometric center before it's rejected as miscentered
CENTRE_TOL = {
    "Core": 3,
    "Intermediate": 999,
    "Outer_1": 999,
    "Outer_2": 999,
}

MIN_VALID_FRACTION = 0.4

# Which patch(es) feed each stack (currently each patch stacks alone)
STACK_GROUPS = {
    "Core": ["Core"],
    "Intermediate": ["Intermediate"],
    "Outer_1": ["Outer_1"],
    "Outer_2": ["Outer_2"],
}

# ─────────────────────────────────────────────────────────────────────────
# Radial profile
#
# STITCH_RADII fixes the radius at which each patch hands off to the
# next one (keyed by the *inner* patch's name — e.g. "Core": 40 means
# "use Core out to r=40px, then switch to Intermediate"). At that
# radius, radial_profile.py interpolates both patches' profiles and
# forces them to match exactly there, so there's no visible jump — no
# manual y-shift constant needed, just this one radius per boundary.
#
# If a patch has no stacked data (0 stars survived stacking.py), it's
# skipped and the boundary connects the two remaining neighbors
# directly, still at the inner one's STITCH_RADII entry.
#
# A boundary with no entry here (e.g. if you add a 5th patch) falls
# back to auto-detecting the switch point as the start of the two
# patches' overlapping valid-radius range — see STITCH_EDGE_MARGIN_FRAC
# / STITCH_MIN_OVERLAP_POINTS in radial_profile.py's _auto_stitch_radius().
# ─────────────────────────────────────────────────────────────────────────
COLOR_MAP = {
    "Core": "teal",
    "Intermediate": "royalblue",
    "Outer_1": "darkorange",
    "Outer_2": "firebrick",
}

STITCH_RADII = {
    "Core": 40,
    "Intermediate": 90,
    "Outer_1": 250,
}

# A patch's profile is only trusted out to this fraction of its cutout
# half-size, to stay clear of Cutout2D edge/corner artifacts. Used both
# as a sanity check on STITCH_RADII and by the overlap-based auto-stitch
# fallback in radial_profile.py.
STITCH_EDGE_MARGIN_FRAC = 0.45

# Minimum number of overlapping radius bins required between two
# patches for the auto-stitch fallback to trust that overlap; below
# this it raises instead of guessing — see _auto_stitch_radius.
STITCH_MIN_OVERLAP_POINTS = 5

RADIAL_PROFILE_MAX = 1500  # px, outer edge of the radial profile


def stacked_psf_path(patch_name: str) -> str:
    """Path to a patch's median-stacked PSF FITS file."""
    return os.path.join(patch_dir(patch_name), "stacked", f"median_psf_{patch_name}.fits")


def sigma_psf_path(patch_name: str) -> str:
    """Path to a patch's per-pixel uncertainty (sigma) stack FITS file, from stacking.py."""
    return os.path.join(patch_dir(patch_name), "stacked", f"sigma_psf_{patch_name}.fits")


if __name__ == "__main__":
    scaffold_project()

    # Local import to avoid a circular import (helpers.py imports this
    # module). Asking here — right after scaffolding — means a user
    # running `python -m polaris.config` on its own gets the same
    # intermediate-plots choice that run_pipeline.py asks up front,
    # instead of only being asked the first time a step tries to plot.
    from .helpers import setup_plot_mode
    setup_plot_mode()