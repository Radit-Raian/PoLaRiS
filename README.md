# PoLaRiS: Point-spread Function Reconstruction of Stellar Sources

PoLaRiS is an automated pipeline for PSF modeling and reconstruction from
multi-band astronomical imaging data, currently optimized for the
Kilo-Degree Survey (KiDS) and the DESI Legacy Imaging Surveys. It takes raw
survey tiles all the way through source detection, Gaia validation, star
selection, cutout extraction, masking, stacking, and radial-profile
construction, producing one continuous, stitched PSF profile per tile set.

## Table of Contents

1. [Step 0: Getting Started](#step-0-getting-started)
2. [Step 1: Configuration](#step-1-configuration)
3. [Step 2: Data Preprocessing](#step-2-data-preprocessing)
4. [Step 3: SExtractor Setup and Run](#step-3-sextractor-setup-and-run)
5. [Step 4: Background Subtraction](#step-4-background-subtraction)
6. [Step 5: Point Source Identification](#step-5-point-source-identification)
7. [Step 6: Gaia Validation](#step-6-gaia-validation)
8. [Step 7: Star Selection](#step-7-star-selection)
9. [Step 8: Cutout Generation](#step-8-cutout-generation)
10. [Step 9: Masking](#step-9-masking)
11. [Step 10: Stacking](#step-10-stacking)
12. [Step 11: Radial Profile Construction](#step-11-radial-profile-construction)
13. [Complete File Reference](#complete-file-reference)
14. [Acknowledgements](#acknowledgements)

---

## Step 0: Getting Started

**Purpose:** Install every dependency the pipeline needs — the Python
package itself and the SExtractor binary it shells out to.

### 0.1 Clone the repository

```bash
git clone https://github.com/Radit-Raian/PoLaRiS.git
cd PoLaRiS
pip install -r requirements.txt
```

### 0.2 Build SExtractor

SExtractor is a compiled C binary, not a Python package, so it is built
and installed separately:

```bash
git clone https://github.com/astromatic/sextractor
cd sextractor
sh autogen.sh
./configure
make -j
sudo make install
sex
```

If it installs somewhere other than `/usr/local/bin/sex`, point
`POLARIS_SEX_BIN` at it in Step 1 rather than editing code.

---

## Step 1: Configuration

**File:** `polaris/config.py`
**Purpose:** Define every path and per-patch parameter used by every other
module, in one place, so a single edit is picked up pipeline-wide.

### Directory layout

Paths are auto-detected relative to `config.py` (`PROJECT_ROOT` is one
level up from the file), so cloning the repo anywhere and dropping data
into `input/` works with zero path edits:

```
<project root>/
├── input/
│   ├── FITS_images/   survey tiles, e.g. KIDS_9.3_-31.2.fits
│   └── sextractor/    gauss_*.conv, default.nnw/.psf/.param
├── output/
│   ├── se_in_out/         SExtractor catalogs, segmaps, backgrounds
│   ├── patches/<name>/    cutouts / masks / stacked, per patch
│   ├── gaia_cache/<tile>/ cached Gaia DR3 query results (auto-fetched)
│   └── pipeline_plots/<step>/   saved figures
└── polaris/           this package
```

Run the scaffolder once to create the empty `input/`/`output/` skeleton:

```bash
python -m polaris.config
```

This also asks whether pipeline plots should display interactively or
save to `output/pipeline_plots/<step>/`:

```
Show intermediate plots interactively? [y/N]:
```

Answering `N` (or running non-interactively) saves every step's figures
to disk instead of popping up a window — useful on a remote server or
when processing many tiles unattended.

---

## Step 2: Data Preprocessing

**File:** `polaris/preprocessing.py`
**Purpose:** Visualize every raw FITS tile before any downstream
processing runs, so data-quality issues are caught early.

`inspect_raw_fits()` downsamples each tile (`downsample=4` by default)
and displays it percentile-scaled (1st–99th) on a dark background, in a
grid with one panel per tile, titled from the FITS header's `OBJECT`
keyword.

```bash
python -m polaris.preprocessing
```

**Output:** `pipeline_plots/preprocessing/raw_fits_preview.png`

---

## Step 3: SExtractor Setup and Run

**File:** `polaris/run_sextractor.py`
**Purpose:** Write a per-tile SExtractor configuration, run source
detection, and preview the resulting segmentation maps.

### Per-tile parameter resolution

Three values are read straight from each tile's FITS header rather than
hardcoded, since they vary tile to tile:

| Header keyword | SExtractor parameter | Fallback |
|---|---|---|
| `PSF_RAD` | `SEEING_FWHM` (arcsec) | `0.7` |
| `STATMAX × 0.95` | `SATUR_LEVEL` (95% of brightest pixel) | `3.5e-8 × 0.95` |
| `CD1_1` → `CDELT1` → `config.DEFAULT_PIXSCALE` | `PIXEL_SCALE` (arcsec/px) | `0.2` |

SExtractor is then run via `subprocess`, and its ASCII catalog is parsed
into a DataFrame to report how many detections passed `CLASS_STAR ≥ 0.84`.

```bash
python -m polaris.run_sextractor
```

### Output files

| File | Contents |
|---|---|
| `se_in_out/outparams_<tile>.cat` | `NUMBER, MAG_AUTO, MAGERR_AUTO, X_IMAGE, Y_IMAGE, FLAGS, ELLIPTICITY, CLASS_STAR, FLUX_RADIUS` |
| `se_in_out/segmap_<tile>.fits` | Segmentation map |
| `se_in_out/bkg_<tile>.fits` | Background map |
| `se_in_out/bkg_rms_<tile>.fits` | Background RMS map |
| `pipeline_plots/run_sextractor/segmaps_preview.png` | Grid preview of every tile's segmentation map |

---

## Step 4: Background Subtraction

**File:** `polaris/background.py`
**Purpose:** Visualize the raw image, SExtractor's background model,
its RMS, and the background-subtracted result, side by side, using the
checkimages Step 3 produced.

`plot_background_subtraction()` builds one row per tile with four
columns — Original, Background, Background RMS, Background Subtracted
(`raw − background`) — each downsampled and percentile-scaled.

```bash
python -m polaris.background
```

**Output:** `pipeline_plots/background/background_subtraction.png`

---

## Step 5: Point Source Identification

**File:** `polaris/point_sources.py`
**Purpose:** Select point-like detections and visualize the
FLUX_RADIUS–MAG_AUTO stellar locus.

`_read_sextractor_catalog()` parses a SExtractor ASCII_HEAD catalog's `#`
header lines to reconstruct column names, so every downstream step reads
the same columns consistently.

`plot_point_sources()` overlays detections with `CLASS_STAR ≥ 0.84` (a
fixed threshold) on each tile's downsampled image.

`flux_radius_vs_mag()` instead derives a **per-tile adaptive** threshold:
it histograms `CLASS_STAR` (100 bins) and finds the valley between `lo=0.3`
and `hi=0.9` — the natural galaxy/star separation point — clipped to
`[0.5, 0.9]`. It then concatenates every tile's catalog and point-source
table and plots FLUX_RADIUS vs. MAG_AUTO for both.

```bash
python -m polaris.point_sources
```

### Output files

| File | Contents |
|---|---|
| `pipeline_plots/point_sources/point_sources.png` | Fixed-threshold (0.84) overlay per tile |
| `pipeline_plots/point_sources/flux_radius_vs_mag.png` | All sources vs. adaptively-selected point sources |

`flux_radius_vs_mag()` returns `(catalog_all, stars_all)` for reuse by
later steps.

---

## Step 6: Gaia Validation

**File:** `polaris/gaia_validation.py`
**Purpose:** Cross-match SExtractor detections against Gaia DR3 and
assign each surviving star its authoritative patch label.

### Pipeline, per tile

1. Select PSF-quality sources: `CLASS_STAR ≥ 0.84`, `FLAGS ≤ 7`.
2. Split off **Outer_2** (the faintest region, by `config.MAG_CUTS`) as
   SExtractor-only — these stars are too faint and numerous to rely on
   Gaia for.
3. Query Gaia DR3 within the tile's sky footprint (`ra`, `dec` bounds
   from the four corner WCS coordinates), restricted to `G ≥ 11` and
   either a high stellar-classification probability or `G < 14`. Results
   are cached to `output/gaia_cache/<tile>/gaia_sources.csv` so repeat
   runs don't re-query.
4. Cross-match via a kd-tree nearest-neighbor search on unit-sphere
   Cartesian coordinates. The match radius is not a fixed arcsec cut —
   it's found per tile from the separation histogram: the first bin
   where the count drops below `max(3, 5% of the peak bin)`.
5. Combine the Gaia-matched and SE-only (Outer_2) sources, drop
   duplicate `NUMBER`s (keeping the smallest separation), then reject
   sources whose brightest pixel isn't within 3 px of the cutout center
   — saturated sources (`FLAGS & 4`) are exempted, since bleed can
   genuinely shift their peak off-center.
6. Label every surviving star's `REGION` via magnitude alone
   (`config.MAG_CUTS`, using Gaia G when available, calibrated
   `MAG_AUTO` otherwise) — **this is the pipeline's authoritative
   per-star patch assignment.**

```bash
python -m polaris.gaia_validation
```

### `master.csv` columns

`NUMBER, MAG_AUTO, MAGERR_AUTO, X_IMAGE, Y_IMAGE, FLAGS, ELLIPTICITY,
CLASS_STAR, FLUX_RADIUS, RA, Dec, MAG_CAL, GAIA_RA, GAIA_Dec, GAIA_Gmag,
Sep_arcsec, TILE, SOURCE, REGION`

### Output files

| File | Contents |
|---|---|
| `output/gaia_validation/master.csv` | Combined catalog across all tiles |
| `output/gaia_validation/<tile>.csv` | Per-tile matched catalog |
| `output/gaia_cache/<tile>/gaia_sources.csv` | Raw cached Gaia query result |
| `pipeline_plots/gaia_validation/source_overlay.png` | Region-coded source overlay per tile |

---

## Step 7: Star Selection

**File:** `polaris/star_selection.py`
**Purpose:** Visually confirm that the flux-radius-based patch
boundaries (`config.PATCHES`) actually capture the intended stellar
locus in `master.csv`.

`plot_star_selection()` scatters `FLUX_RADIUS` vs. `MAG_AUTO` for every
matched star, then overlays each patch's rectangle from `config.PATCHES`
with its star count in the legend — the split Intermediate range is
drawn as one labeled box plus one unlabeled box so the legend isn't
duplicated.

```bash
python -m polaris.star_selection
```

**Output:** `pipeline_plots/star_selection/star_selection.png`

---

## Step 8: Cutout Generation

**File:** `polaris/cutouts.py`
**Purpose:** Extract a fixed-size cutout around each selected star,
sized by which patch it falls into.

`extract_cutout()` is the single source of truth for turning
`(tile image, x, y, size)` into a saved array, using `Cutout2D` with
`mode='partial'` so stars near a tile edge still produce a
correctly-shaped array (filled outside the tile). `masking.py` imports
this same function directly, so the two steps' extraction logic can
never drift apart.

`generate_cutouts()` loops over every matched star, determines its
patch via `helpers.get_patch()` against `config.PATCHES`, rejects stars
too close to the tile edge (`edge_margin_frac=0.55` of the cutout size),
extracts and validates the cutout (`is_valid_cutout()` rejects >20% NaN/
zero pixels, `config.MAX_EMPTY_FRACTION`), and saves it as
`<tile>_<star_number>_<size>.fits` under
`output/patches/<folder>/cutouts/`.

```bash
python -m polaris.cutouts
```

### Output files

| File | Contents |
|---|---|
| `patches/<folder>/cutouts/<tile>_<number>_<size>.fits` | Per-star science cutout |
| `pipeline_plots/cutouts/cutout_preview_<folder>.png` | 5 random cutouts per patch |

---

## Step 9: Masking

**File:** `polaris/masking.py`
**Purpose:** Build a contamination mask for every star and produce all
five per-star products needed for stacking.

This step wipes and regenerates its own cutouts from scratch (it does
not reuse Step 8's output on disk) because the segmap cutout must be
extracted with the exact same center and size as the science cutout,
and centering can shift once the target-object check below is applied.

For every star:

1. Extract the science cutout and the matching segmap cutout (via
   `cutouts.extract_cutout`).
2. Identify the target object's segmentation ID at the cutout's center
   pixel; skip the star if the center isn't on any detected object.
3. Build a **primary mask**: star + background pixels = 0, everything
   else = 1.
4. Flag every segmented pixel belonging to a *different* object, then
   grow that flag outward with binary dilation (`DILATION_ITERATIONS=4`,
   full 8-connectivity) so faint wings of nearby contaminants are also
   excluded, not just their detected core.
5. Build the **final mask** as the inverse of the primary mask, and the
   **multiplied** image as `cutout × final_mask`.

```bash
python -m polaris.masking
```

### Output files (per star, under `patches/<folder>/`)

| Subfolder | Contents |
|---|---|
| `cutouts/` | Science cutout |
| `segmaps/` | Matching segmentation-map cutout |
| `primary_masking/` | Star+background=0, contaminants=1 |
| `final_masking/` | Inverse of primary — 1=keep, 0=mask out |
| `multiplied/` | `cutout × final_mask` |

Plus `pipeline_plots/masking/mask_preview_<folder>_<star_id>.png` — a
5-panel preview (Cutout, Segmap, Primary Mask, Final Mask, Masked Image)
for 5 random stars per patch.

---

## Step 10: Stacking

**File:** `polaris/stacking.py`
**Purpose:** Reject bad cutouts, decorrelate non-circular PSF features,
normalize, and stack every patch's surviving stars.

### Rejection criteria, applied in sequence

| Check | Threshold | Rejects |
|---|---|---|
| Shape | `!= expected_size` | Malformed cutouts |
| Valid pixel fraction | `< MIN_VALID_FRACTION` (0.4) | Too heavily masked/edge-affected |
| Masked fraction | `> config.MAX_MASKED_FRAC[patch]` | Same check, per-patch tolerance |
| Peak offset | `> config.CENTRE_TOL[patch]` px | Miscentered stars |
| Annulus normalization | non-finite or ≤ 0 | Bad flux normalization |

| Patch | Norm. annulus (px) | Max masked fraction | Centre tolerance (px) |
|---|---|---|---|
| Core | 10 – 15 | 0.20 | 3 |
| Intermediate | 100 – 140 | 0.30 | 999 (unconstrained) |
| Outer_1 | 200 – 220 | 0.50 | 999 (unconstrained) |
| Outer_2 | 400 – 410 | 0.80 | 999 (unconstrained) |

Surviving cutouts are rotated by a **random continuous angle** in
`[0, 360)` (bilinear interpolation, `cval=NaN`). A continuous angle is
required — 90°-only rotation would leave pixel-grid-aligned features
like diffraction spikes and chip seams uncorrelated across stars.

Each rotated cutout is normalized by its annulus median flux, then
stacked into a median stack, a mean stack, and a per-pixel uncertainty
map (Garate et al.): `(P84 − P16) / 2 × 1.253 / √n_contributing`.

```bash
python -m polaris.stacking
```

A patch with zero accepted stars is skipped entirely — no FITS files
are written for it, and it's logged as `[SKIP] <patch>: no valid PSFs`.

### Output files

| File | Contents |
|---|---|
| `patches/<patch>/stacked/median_psf_<patch>.fits` | Median stack |
| `patches/<patch>/stacked/mean_psf_<patch>.fits` | Mean stack |
| `patches/<patch>/stacked/sigma_psf_<patch>.fits` | Per-pixel uncertainty |
| `pipeline_plots/stacking/stack_<patch>.png` | Log-stretched median/mean side by side, with normalization-annulus rings drawn |

A summary table (Total / Accepted / rejection counts per category) is
printed for every patch after the run.

---

## Step 11: Radial Profile Construction

**File:** `polaris/radial_profile.py`
**Purpose:** Compute each patch's radial profile from its median stack,
then stitch them into one continuous curve — automatically.

`available_patches()` checks which patches actually have a median stack
on disk; a patch missing one (0 accepted stars in Step 10) is printed as
`[SKIP]` and excluded, so the profile is built from whichever patches
exist even if that's just one.

`compute_patch_profiles()` radial-profiles every available patch over
one shared radius grid (linear from 0–10 px, then log-spaced from 10 px
to `config.RADIAL_PROFILE_MAX`), normalized to peak = 1.

### Auto-stitching — no manual y-shift

For each pair of adjacent patches (ordered by cutout size):

1. A radius is **valid** for a patch if its profile is finite, positive,
   and inside `config.STITCH_EDGE_MARGIN_FRAC` of its cutout — clear of
   `Cutout2D` edge/corner artifacts. No SNR comparison is involved.
2. The stitch radius comes from `config.STITCH_RADII` when the inner
   patch has an entry there (manual); otherwise it's auto-detected as
   the **start** of the two patches' overlapping valid-radius range —
   the earliest point the outer patch can be trusted, so the inner
   patch's better-resolved core is kept as long as possible.
3. Both profiles are interpolated (log-log) at that exact radius, and
   the outer patch is scaled so the two lines match there exactly. That
   scale factor is the only "shift" applied — it's derived from the
   data, not chosen by hand.
4. If a boundary has fewer than `config.STITCH_MIN_OVERLAP_POINTS`
   overlapping radii and no manual entry, it raises rather than
   guessing — add a `config.STITCH_RADII` entry for that boundary.

```bash
python -m polaris.radial_profile
```

`print_stitch_report()` prints the derived scale, stitch radius, source
("manual"/"auto-overlap"), overlap size, and a `clamped` flag per
boundary, so the automatic choice can be sanity-checked before trusting
it.

**Output:** `pipeline_plots/radial_profile/radial_profile.png` — the
continuous PSF profile (log-log), with each patch's contributing radius
range shaded.

---

## Complete File Reference

| Step | Module | Key outputs |
|---|---|---|
| 1 | `config.py` | `input/`, `output/` folder skeleton |
| 2 | `preprocessing.py` | `raw_fits_preview.png` |
| 3 | `run_sextractor.py` | `outparams_<tile>.cat`, `segmap_<tile>.fits`, `bkg_<tile>.fits`, `bkg_rms_<tile>.fits` |
| 4 | `background.py` | `background_subtraction.png` |
| 5 | `point_sources.py` | `point_sources.png`, `flux_radius_vs_mag.png` |
| 6 | `gaia_validation.py` | `master.csv`, per-tile CSVs, `source_overlay.png` |
| 7 | `star_selection.py` | `star_selection.png` |
| 8 | `cutouts.py` | `patches/<folder>/cutouts/*.fits`, `cutout_preview_<folder>.png` |
| 9 | `masking.py` | `patches/<folder>/{cutouts,segmaps,primary_masking,final_masking,multiplied}/*.fits` |
| 10 | `stacking.py` | `patches/<patch>/stacked/{median,mean,sigma}_psf_<patch>.fits` |
| 11 | `radial_profile.py` | `radial_profile.png` |

---

## Acknowledgements

- SExtractor: Bertin, E. & Arnouts, S. (1996), *Astronomy and Astrophysics
  Supplement Series*, 117, 393–404
- KiDS Survey: <https://kids.strw.leidenuniv.nl/>
- DESI Legacy Surveys: <https://www.legacysurvey.org/>

### References

| # | Title | Link |
|---|---|---|
| 1 | The Buildup of the Intracluster Light of A85 as Seen by Subaru's Hyper Suprime-Cam | IOP Science |
| 2 | The Hyper Suprime-Cam extended point spread functions and applications | MNRAS |

---

*Pipeline documentation for PSF Modeling*
