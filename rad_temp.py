"""
One-off: build the continuous PSF radial profile using only the three
patches that actually produced a stack this run (Core, Intermediate,
Outer_1) — Outer_2 had 0 accepted stars on this single-tile test run.

Run from the repo root:  python radial_profile_3patch.py
"""

import numpy as np

from polaris import config
from polaris.radial_profile import compute_patch_profiles, plot_continuous_profile

PATCHES = ["Core", "Intermediate", "Outer_1"]  # Outer_2 skipped — no stack produced

stack_paths = {name: config.stacked_psf_path(name) for name in PATCHES}
profiles = compute_patch_profiles(stack_paths)

r_core, p_core = profiles["Core"]["r"], profiles["Core"]["profile"] * config.YSHIFT["Core"]
r_int, p_int = profiles["Intermediate"]["r"], profiles["Intermediate"]["profile"] * config.YSHIFT["Intermediate"]
r_o1, p_o1 = profiles["Outer_1"]["r"], profiles["Outer_1"]["profile"] * config.YSHIFT["Outer_1"]

mask_core = r_core <= config.CORE_MAX
mask_int = (r_int > config.CORE_MAX) & (r_int <= config.INT_MAX)
mask_o1 = r_o1 > config.INT_MAX  # goes all the way out since there's no Outer_2 to hand off to

r_plot = np.concatenate([r_core[mask_core], r_int[mask_int], r_o1[mask_o1]])
p_plot = np.concatenate([p_core[mask_core], p_int[mask_int], p_o1[mask_o1]])

valid = np.isfinite(r_plot) & np.isfinite(p_plot) & (r_plot > 0) & (p_plot > 0)
r_plot, p_plot = r_plot[valid], p_plot[valid]

print(f"Stitched profile: {len(r_plot)} points, r={r_plot.min():.1f}-{r_plot.max():.1f} px")

plot_continuous_profile(r_plot, p_plot)