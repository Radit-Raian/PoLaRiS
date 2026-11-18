"""
Run the full PoLaRiS pipeline end to end, in order:

  preprocessing -> run_sextractor -> background -> point_sources
  -> gaia_validation -> star_selection -> cutouts -> masking
  -> stacking -> radial_profile

Each step can also be run independently, e.g.:

    python -m polaris.stacking
    python -m polaris.radial_profile

which is usually more useful while iterating on a single step, since
this script re-runs everything from scratch.

Usage:
    python run_pipeline.py                  # run every step
    python run_pipeline.py --from stacking   # resume from a given step
    python run_pipeline.py --skip-plots      # suppress interactive plots
"""

import argparse

from polaris import (
    preprocessing,
    run_sextractor,
    background,
    point_sources,
    gaia_validation,
    star_selection,
    cutouts,
    masking,
    stacking,
    radial_profile,
)

STEPS = [
    ("preprocessing", preprocessing.main),
    ("run_sextractor", run_sextractor.main),
    ("background", background.main),
    ("point_sources", point_sources.main),
    ("gaia_validation", gaia_validation.main),
    ("star_selection", star_selection.main),
    ("cutouts", cutouts.main),
    ("masking", masking.main),
    ("stacking", stacking.main),
    ("radial_profile", radial_profile.main),
]


def main():
    parser = argparse.ArgumentParser(description="Run the PoLaRiS PSF pipeline.")
    parser.add_argument("--from", dest="start_from", default=None,
                         choices=[name for name, _ in STEPS],
                         help="Resume the pipeline from this step (skips everything before it).")
    args = parser.parse_args()

    start_idx = 0
    if args.start_from:
        start_idx = next(i for i, (name, _) in enumerate(STEPS) if name == args.start_from)

    for name, step_main in STEPS[start_idx:]:
        print(f"\n{'=' * 70}\nSTEP: {name}\n{'=' * 70}")
        step_main()


if __name__ == "__main__":
    main()
