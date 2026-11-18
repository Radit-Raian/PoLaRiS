"""
Step: Star Selection.

Loads the Gaia-matched master catalog and overlays the patch
definitions (Core / Intermediate / Outer_1 / Outer_2) on a
FLUX_RADIUS vs MAG_AUTO plot, to visually confirm the patch
boundaries in config.PATCHES capture the intended stellar locus.

TODO (unresolved): gaia_validation.py now writes a REGION column to
master.csv, assigned purely by magnitude (config.MAG_CUTS) — that's
the pipeline's actual, authoritative per-star patch assignment.
This module's flux-radius-based `_PLOT_REGIONS` (and config.PATCHES,
which cutouts.py currently uses) is a separate, independently-defined
scheme that will disagree with REGION for some stars. Before running
cutouts.py on real data, decide which assignment is correct and make
cutouts.py use master.csv's REGION column directly instead of
re-deriving a patch from flux_radius/mag — see config.py for the same
flag.
"""

from matplotlib import patches as mpatches
import matplotlib.pyplot as plt
import pandas as pd

from . import config
from .helpers import finalize_plot, setup_plot_mode

# (x_min, x_max, y_min, y_max, color, label) — visual regions matching
# config.PATCHES, with the split Intermediate range drawn as one labeled
# box plus one unlabeled box so the legend isn't duplicated.
_PLOT_REGIONS = [
    (0, 3.25, 16, 18, "teal", "Core"),
    (0, 4, 14.5, 15, "royalblue", "Intermediate"),
    (2, 8, 12.5, 14.5, "royalblue", None),
    (8, 25, 11, 13.25, "darkorange", "Outer_1"),
    (25, 175, 8, 11, "firebrick", "Outer_2"),
]


def plot_star_selection(matched_catalog_path=None):
    """Scatter FLUX_RADIUS vs MAG_AUTO for all matched stars, with patch regions overlaid."""
    matched_catalog_path = matched_catalog_path or config.MATCHED_CATALOG_PATH
    matched_catalog = pd.read_csv(matched_catalog_path)

    fig, ax = plt.subplots(figsize=(10, 7))

    ax.scatter(matched_catalog["FLUX_RADIUS"], matched_catalog["MAG_AUTO"], s=20, color="gray", alpha=0.5)

    for x_min, x_max, y_min, y_max, color, label in _PLOT_REGIONS:
        rect = mpatches.Rectangle(
            (x_min, y_min), x_max - x_min, y_max - y_min,
            linewidth=0, edgecolor=None, facecolor=color, alpha=0.12,
        )
        ax.add_patch(rect)

        mask = (
            (matched_catalog["FLUX_RADIUS"] >= x_min) & (matched_catalog["FLUX_RADIUS"] <= x_max)
            & (matched_catalog["MAG_AUTO"] >= y_min) & (matched_catalog["MAG_AUTO"] <= y_max)
        )
        count = mask.sum()

        ax.scatter(
            matched_catalog.loc[mask, "FLUX_RADIUS"], matched_catalog.loc[mask, "MAG_AUTO"],
            s=35, color=color, alpha=0.9, label=f"{label}: {count}" if label else "_nolegend_",
        )

    ax.set_xlim(-3, 25)
    ax.set_ylim(matched_catalog["MAG_AUTO"].max() + 1, matched_catalog["MAG_AUTO"].min() - 1)
    ax.set_xlabel("FLUX_RADIUS (pixels)")
    ax.set_ylabel("MAG_AUTO")
    ax.set_title("FLUX_RADIUS vs MAG_AUTO")
    ax.legend()
    ax.grid(True)
    ax.set_axisbelow(True)

    plt.tight_layout()
    finalize_plot(fig, "star_selection.png", subdir="star_selection")


def main():
    plot_star_selection()


if __name__ == "__main__":
    main()