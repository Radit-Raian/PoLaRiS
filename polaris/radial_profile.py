"""
Step: Radial Profile Construction.

Computes each patch's radial profile from its median-stacked PSF FITS
file, then stitches the patch profiles into one continuous curve.

Stitching is automatic. For each pair of adjacent patches (ordered by
cutout size — Core, Intermediate, Outer_1, Outer_2):

  1. Both patches' flux profiles are radial-profiled over the *same*
     full radius grid (from stacking.py's median stacks).
  2. A radius is "valid" for a patch if its profile is finite/positive
     and comfortably inside its cutout (config.STITCH_EDGE_MARGIN_FRAC,
     to stay clear of Cutout2D edge/corner artifacts) — no SNR
     involved.
  3. The overlap between two patches' valid radii gives:
       - a stitch radius: the *start* of that overlap — the first
         (innermost) radius where both patches' profiles are valid —
         found automatically from the two profiles themselves, no
         hand-picked constant.
       - a scale factor: both profiles are interpolated at that exact
         stitch radius and the outer patch is scaled so the two lines
         match there exactly — this is the only "y-shift", and it's
         derived, not chosen.
  4. If a boundary has too little overlap (config.
     STITCH_MIN_OVERLAP_POINTS), it raises rather than guessing —
     add a config.STITCH_RADII entry for that boundary instead.

No manual y-shift constant is needed, and no SNR/sigma-stack
comparison is involved; changing star counts, seeing, or magnitude
cuts re-derives the stitch radius and scale automatically on the next
run. See auto_stitch_profiles() for the implementation and
STITCH_REPORT below for how to sanity-check what it picked.

Not every patch is guaranteed to have a stack — stacking.py skips a
patch entirely (writing no FITS files) if zero stars survive its
rejection cuts, which happens on tiles too small/shallow to have
enough very bright stars for e.g. Outer_2. available_patches() checks
for that up front, so a missing patch is skipped with a clear message
instead of crashing, and the continuous profile is built from
whichever patches actually exist — even if that's just one.
"""

import os

import numpy as np
from astropy.io import fits
from photutils.profiles import RadialProfile
import matplotlib.pyplot as plt

from . import config
from .helpers import finalize_plot, setup_plot_mode


def available_patches(patch_order=None):
    """
    Return the subset of `patch_order` that actually has a median
    stack on disk, in the same relative order.

    A patch is missing its stack when stacking.py accepted zero stars
    for it (it prints "[SKIP] <patch>: no valid PSFs" and writes
    nothing) — not a bug, just a patch with no usable data on this
    tile. Prints one line per skipped patch so it's obvious why the
    final profile is missing a region, rather than surfacing as a
    FileNotFoundError.
    """
    patch_order = patch_order or list(config.PATCH_SIZES.keys())

    available = []
    for name in patch_order:
        med_path = config.stacked_psf_path(name)
        if not os.path.exists(med_path):
            print(f"[SKIP] {name}: missing median stack — "
                  f"likely 0 accepted stars in stacking.py. Excluded from radial profile.")
        else:
            available.append(name)

    return available


def compute_patch_profiles(patch_order=None, stack_paths=None):
    """
    Compute each patch's radial profile, over one shared radius grid
    so patches can be directly compared/overlapped.

    Flux profile is normalized to peak=1 (method="max").

    Patches missing a stack (see available_patches()) are silently
    excluded rather than raising. Raises only if none are available.

    Returns {patch_name: {"r", "profile"}}.
    """
    patch_order = available_patches(patch_order)
    if not patch_order:
        raise RuntimeError(
            "No patches have a stacked PSF yet — run stacking.py first, "
            "and check its per-patch summary table for 0-accepted patches."
        )

    stack_paths = stack_paths or {name: config.stacked_psf_path(name) for name in patch_order}

    edge_radii = np.concatenate([
        np.arange(0, 10, 1),
        np.logspace(np.log10(10), np.log10(config.RADIAL_PROFILE_MAX), 300)[1:],
    ])

    profiles = {}
    for name in patch_order:
        psf = fits.getdata(stack_paths[name]).astype("float64")
        ny, nx = psf.shape

        rp_flux = RadialProfile(psf, (nx / 2, ny / 2), edge_radii)
        raw_profile = rp_flux.profile.copy()
        peak = np.nanmax(raw_profile)
        flux_profile = raw_profile / peak

        profiles[name] = {"r": rp_flux.radius, "profile": flux_profile}

    return profiles


def _valid_mask(entry, cutout_size, r, edge_margin_frac):
    """Radii where a patch's profile is usable: finite, positive, and inside its cutout."""
    size_limit = cutout_size * edge_margin_frac
    finite_positive = np.isfinite(entry["profile"]) & (entry["profile"] > 0)
    return finite_positive & (r <= size_limit)


def _interp_at_radius(r, profile, r_target):
    """
    Interpolate `profile` (as a function of `r`) at `r_target`, in
    log-log space (profiles are roughly power-law in r, so this is far
    more accurate than linear interpolation near a steep core).

    Only finite, positive points are used. If `r_target` falls outside
    the covered range, np.interp clamps to the nearest edge value
    rather than extrapolating — the caller is warned when that happens.
    """
    mask = np.isfinite(profile) & (profile > 0)
    r_valid, p_valid = r[mask], profile[mask]
    clamped = not (r_valid.min() <= r_target <= r_valid.max())
    value = float(np.exp(np.interp(np.log(r_target), np.log(r_valid), np.log(p_valid))))
    return value, clamped


def _auto_stitch_radius(prev_entry, entry, prev_size, cur_size, r, edge_margin_frac, min_overlap_points):
    """
    Fallback for a boundary with no config.STITCH_RADII entry.

    Finds the overlap between the two patches' valid-radius ranges
    (see _valid_mask — finite/positive, inside the edge margin, no SNR
    involved) and uses its *start* — the first, innermost radius where
    both patches are valid — as the stitch radius. That's the earliest
    point the outer (bigger-cutout) patch's profile can be trusted, so
    handing off there keeps as much of the inner patch's better-resolved
    core as possible.
    """
    overlap = (
        _valid_mask(prev_entry, prev_size, r, edge_margin_frac)
        & _valid_mask(entry, cur_size, r, edge_margin_frac)
    )
    if overlap.sum() < min_overlap_points:
        raise ValueError(
            "No usable overlap and no config.STITCH_RADII entry for this boundary — "
            "check that both patches' stacks cover a shared radius range, or add a manual radius."
        )

    overlap_r = r[overlap]
    r_stitch = float(overlap_r[0])
    return r_stitch, int(overlap.sum())


def auto_stitch_profiles(profiles, patch_order=None, patch_sizes=None, stitch_radii=None,
                          edge_margin_frac=None, min_overlap_points=None):
    """
    Stitch adjacent patches into one continuous (r, profile) curve.

    The switch radius between two patches comes from `stitch_radii`
    (config.STITCH_RADII by default, keyed by the inner patch's name)
    when available; otherwise it's auto-detected as the start of the
    two patches' overlapping valid-radius range (see
    _auto_stitch_radius) — no SNR/sigma-stack comparison involved.

    Either way, the scale factor is derived by interpolating BOTH
    patches' profiles at that exact switch radius and forcing them to
    match there — "match the lines" at the boundary, no manual y-shift.

    Returns (r_plot, p_plot, report) where `report` is a list of dicts
    (one per boundary) with the inner/outer patch names, the derived
    scale, the stitch radius, its source ("manual"/"auto-overlap"), the
    overlap size behind an auto radius, and a `clamped` flag if the
    radius fell outside a patch's covered range — print via
    print_stitch_report() to sanity-check before trusting it.
    """
    patch_order = patch_order or list(config.PATCH_SIZES.keys())
    patch_sizes = patch_sizes or config.PATCH_SIZES
    stitch_radii = stitch_radii if stitch_radii is not None else config.STITCH_RADII
    edge_margin_frac = edge_margin_frac if edge_margin_frac is not None else config.STITCH_EDGE_MARGIN_FRAC
    min_overlap_points = min_overlap_points if min_overlap_points is not None else config.STITCH_MIN_OVERLAP_POINTS

    r = profiles[patch_order[0]]["r"]  # shared grid — every patch used the same edge_radii
    scaled = {patch_order[0]: profiles[patch_order[0]]["profile"].copy()}
    boundaries = []
    report = []

    for prev_name, name in zip(patch_order, patch_order[1:]):
        prev_entry, entry = profiles[prev_name], profiles[name]

        if prev_name in stitch_radii:
            r_stitch, source, n_overlap = stitch_radii[prev_name], "manual", None
        else:
            r_stitch, n_overlap = _auto_stitch_radius(
                prev_entry, entry, patch_sizes[prev_name], patch_sizes[name],
                r, edge_margin_frac, min_overlap_points,
            )
            source = "auto-overlap"

        # Match the lines: interpolate both profiles at r_stitch and
        # solve for the scale that makes the outer patch agree with the
        # (already-scaled) inner one there.
        prev_value, prev_clamped = _interp_at_radius(r, scaled[prev_name], r_stitch)
        cur_value, cur_clamped = _interp_at_radius(r, entry["profile"], r_stitch)
        scale = prev_value / cur_value
        scaled[name] = entry["profile"] * scale

        boundaries.append(r_stitch)
        report.append({
            "inner": prev_name, "outer": name, "scale": scale, "r_stitch": r_stitch,
            "source": source, "n_overlap": n_overlap, "clamped": prev_clamped or cur_clamped,
        })

    # Combine: each patch contributes only within its stitch-derived range.
    lows = [-np.inf] + boundaries
    highs = boundaries + [config.RADIAL_PROFILE_MAX]

    r_pieces, p_pieces = [], []
    for name, lo, hi in zip(patch_order, lows, highs):
        mask = (r > lo) & (r <= hi)
        r_pieces.append(r[mask])
        p_pieces.append(scaled[name][mask])

    r_plot, p_plot = np.concatenate(r_pieces), np.concatenate(p_pieces)
    valid = np.isfinite(r_plot) & np.isfinite(p_plot) & (r_plot > 0) & (p_plot > 0)
    return r_plot[valid], p_plot[valid], report


def print_stitch_report(report):
    """Print the derived scale factor and stitch radius for each boundary, for a quick sanity check."""
    print("\nRadial-profile stitching (auto-derived):")
    for row in report:
        overlap = f"n_overlap={row['n_overlap']}" if row["n_overlap"] is not None else "n_overlap=n/a (manual)"
        flag = "  [clamped — r_stitch outside a patch's covered range]" if row["clamped"] else ""
        print(
            f"  {row['inner']:<13} -> {row['outer']:<13} "
            f"scale={row['scale']:.4e}  r_stitch={row['r_stitch']:7.2f}px  "
            f"source={row['source']:<12} {overlap}{flag}"
        )


def plot_continuous_profile(r_plot, p_plot, report, patch_order=None, color_map=None):
    """Plot the stitched continuous PSF radial profile, with patch regions shaded at their auto-derived boundaries."""
    patch_order = patch_order or list(config.PATCH_SIZES.keys())
    color_map = color_map or config.COLOR_MAP
    boundaries = [row["r_stitch"] for row in report]

    fig = plt.figure(figsize=(8, 6))

    edges = [0] + boundaries + [config.RADIAL_PROFILE_MAX]
    for name, lo, hi in zip(patch_order, edges, edges[1:]):
        plt.axvspan(lo, hi, alpha=0.12, color=color_map[name], label=name)

    plt.plot(r_plot, p_plot, lw=2, color="black")

    plt.xscale("log")
    plt.yscale("log")
    plt.xlim(0.5, config.RADIAL_PROFILE_MAX)
    plt.ylim(1e-11, 2)
    plt.xlabel("Radius (pixels)", fontsize=12)
    plt.ylabel("Normalised Surface Brightness", fontsize=12)
    plt.title("Continuous PSF Radial Profile (auto-stitched)", fontsize=13)
    plt.legend(frameon=False, fontsize=10)
    plt.minorticks_on()
    plt.tight_layout()
    finalize_plot(fig, "radial_profile.png", subdir="radial_profile")


def build_continuous_profile():
    """Full step: compute per-patch profiles for whichever patches exist, auto-stitch, report, and plot."""
    profiles = compute_patch_profiles()
    patch_order = list(profiles.keys())  # only the patches that actually had stacks
    print(f"\nBuilding radial profile from: {', '.join(patch_order)}")

    r_plot, p_plot, report = auto_stitch_profiles(profiles, patch_order=patch_order)
    print_stitch_report(report)
    plot_continuous_profile(r_plot, p_plot, report, patch_order=patch_order)
    return r_plot, p_plot, report


def main():
    build_continuous_profile()


if __name__ == "__main__":
    main()