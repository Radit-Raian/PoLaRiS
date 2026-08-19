# PoLaRiS: Point-spread Function Reconstruction of Stellar Sources

PoLaRiS (Point-spread Function Reconstruction of Stellar Sources) is a modular Python pipeline for constructing extended point-spread function (PSF) profiles from stellar sources in astronomical imaging data. It automates the sequence from source detection through Gaia-based stellar validation, multi-scale cutout generation, contamination masking, PSF stacking, and construction of a continuous radial PSF profile.

The pipeline is currently developed and tested for wide-field optical imaging, with particular application to data from the **Kilo-Degree Survey (KiDS)** and the **DESI Legacy Imaging Surveys**. Its primary goal is to provide a reproducible workflow for reconstructing PSF profiles over a substantially larger radial range than can generally be obtained from a single stellar cutout.

## Table of Contents

- [Scientific Motivation](#scientific-motivation)
- [Repository Structure](#repository-structure)
- [Installation](#installation)
- [Dependencies](#dependencies)
- [Configuration](#configuration)
- [Pipeline](#pipeline)
  1. [FITS Inspection and Preprocessing](#1-fits-inspection-and-preprocessing)
  2. [SExtractor Source Detection](#2-sextractor-source-detection)
  3. [Background Subtraction](#3-background-subtraction)
  4. [Point-source Identification](#4-point-source-identification)
  5. [Gaia Validation](#5-gaia-validation)
  6. [Stellar Selection and Patch Definition](#6-stellar-selection-and-patch-definition)
  7. [Cutout Generation](#7-cutout-generation)
  8. [Contamination Masking](#8-contamination-masking)
  9. [PSF Stacking](#9-psf-stacking)
  10. [Radial Profile Construction](#10-radial-profile-construction)
- [Helper Functions](#helper-functions)
- [Pipeline Execution](#pipeline-execution)
- [Example Notebooks](#example-notebooks)
- [Output Products](#output-products)
- [Reproducibility and Configuration](#reproducibility-and-configuration)
- [Testing](#testing)
- [Current Scope and Limitations](#current-scope-and-limitations)
- [Reference](#reference)
- [Acknowledgements](#acknowledgements)
- [License](#license)
- [Citation](#citation)

## Scientific Motivation

The point-spread function describes the response of an imaging system to an unresolved source. Accurate characterization of the PSF matters for a number of astronomical applications, including the analysis of low-surface-brightness structures and extended diffuse emission.

A PSF measured only from the central region of a stellar image may not adequately describe extended wings at large radii. Constructing an extended PSF therefore requires both suitable stellar sources and measurements over multiple spatial scales. PoLaRiS addresses this by combining stellar sources with different cutout sizes into a single multi-scale PSF reconstruction workflow: detecting sources, identifying and Gaia-validating candidate stars, grouping them by measured properties, extracting cutouts at multiple spatial scales, removing neighbor contamination with segmentation maps, normalizing and stacking the resulting images, and stitching the resulting radial profiles into one continuous curve.

## Repository Structure

PoLaRiS is organized as an installable Python package rather than a sequence of numbered scripts, so each stage is importable, testable, and documented independently:

```
polaris/
├── polaris/
│   ├── __init__.py
│   ├── config.py
│   ├── helpers.py
│   ├── preprocessing.py
│   ├── run_sextractor.py
│   ├── background.py
│   ├── point_sources.py
│   ├── gaia_validation.py
│   ├── star_selection.py
│   ├── cutouts.py
│   ├── masking.py
│   ├── stacking.py
│   └── radial_profile.py
├── run_pipeline.py
├── pyproject.toml
├── README.md
├── LICENSE
├── tests/
└── notebooks/
```

| Module | Description |
|---|---|
| `config.py` | Central configuration of paths, tiles, patch definitions, cutout sizes, normalization parameters, rejection criteria, plotting parameters, and profile-stitching boundaries |
| `helpers.py` | Shared utility functions used throughout the pipeline |
| `preprocessing.py` | FITS-image inspection and downsampling |
| `run_sextractor.py` | SExtractor configuration, source detection, and segmentation-map generation |
| `background.py` | Background subtraction and visualization |
| `point_sources.py` | Point-source identification and adaptive `CLASS_STAR` selection |
| `gaia_validation.py` | Gaia querying and cross-matching |
| `star_selection.py` | Stellar selection and definition of PSF patches |
| `cutouts.py` | Stellar cutout extraction and validation (exposes `generate_cutout()`) |
| `masking.py` | Segmentation-based contamination masking (imports `cutouts.generate_cutout()`) |
| `stacking.py` | Cutout normalization, quality rejection, stacking, and PSF FITS output |
| `radial_profile.py` | Radial-profile construction and multi-patch profile stitching |
| `run_pipeline.py` | Pipeline orchestration and execution entry point |

`cutouts.py` and `masking.py` are kept as separate modules with no duplicated I/O: `cutouts.py` owns "given a star, extract and save a cutout," while `masking.py` owns "given a saved cutout, build masks," calling into `cutouts.py` for the extraction step rather than reimplementing the `Cutout2D` call.

## Installation

### Clone the repository

```bash
git clone https://github.com/Radit-Raian/PoLaRiS.git
cd PoLaRiS
```

### Install the package

```bash
pip install .
```

For development:

```bash
pip install -e .
```

Project dependencies are specified in `pyproject.toml`.

## Dependencies

PoLaRiS uses the scientific Python ecosystem for astronomical image processing and catalog analysis, including:

- NumPy
- SciPy
- Pandas
- Matplotlib
- Astropy
- Photutils
- Astroquery

Exact dependency versions should be obtained from `pyproject.toml`.

### SExtractor

PoLaRiS also requires **SExtractor** for source detection and segmentation-map generation. A source installation may follow:

```bash
git clone https://github.com/astromatic/sextractor
cd sextractor
sh autogen.sh
./configure
make -j
sudo make install
sex
```

The installation procedure may vary depending on the operating system and available system packages.

## Configuration

Pipeline-wide configuration is centralized in `polaris/config.py`, controlling:

- Input and output paths
- `TILE_NAMES`
- `PATCHES` and `PATCH_SIZES`
- `NORM_ANNULUS`
- `MAX_MASKED_FRAC`
- `CENTRE_TOL`
- `COLOR_MAP`
- `YSHIFT`
- `CORE_MAX`, `INT_MAX`, `O1_MAX`

Before processing a new dataset, users should review and adjust these parameters — particularly pixel scale, patch boundaries, cutout sizes, normalization annuli, masking thresholds, centering tolerances, and radial-profile stitching boundaries — since they may depend on the survey, detector sampling, and image characteristics.

## Pipeline

### 1. FITS Inspection and Preprocessing

Implemented in `preprocessing.py`. This stage:

- Loads the input FITS image and inspects the FITS header
- Visualizes the image and downsamples it for rapid inspection
- Sets the pixel scale from the survey

Downsampling is used primarily for visualization and quality inspection of large survey tiles, where displaying a full-resolution image may require substantial memory.

### 2. SExtractor Source Detection

Implemented in `run_sextractor.py`. This stage:

- Defines input/output directories and references the necessary configuration files (`.nnw`, `.psf`, `.param`, convolution filter, etc.)
- Writes a per-tile SExtractor configuration file and loops through all tiles
- Runs SExtractor to detect primary sources and produce source catalogs and segmentation maps
- Visualizes the resulting downsampled segmentation map for validation (kept downsampled to avoid memory issues, since KiDS FITS files are large)

The segmentation maps are later used to identify and remove contamination from neighboring sources.

### 3. Background Subtraction

Implemented in `background.py`. For each tile, the pipeline loads the raw image, the SExtractor background map, and the background RMS map, and computes the background-subtracted image as `I_sub = I_raw − I_background` using Photutils. The stage also provides visualization of the resulting image for inspection.

### 4. Point-source Identification

Implemented in `point_sources.py`. This stage:

- Parses the SExtractor `.cat` header lines to reconstruct column names and loads the catalog into a DataFrame
- Filters for point sources using an initial criterion of `CLASS_STAR ≥ 0.84`, with support for an adaptive `CLASS_STAR` threshold based on the distribution of detected sources
- Overlays candidate sources on the imaging data for visual validation
- Distributes selected stars into four groups by position on the magnitude-vs-flux-radius plot: Core, Intermediate, Outer 1, and Outer 2
- Clips selected stars to acceptable flux radius and magnitude ranges, and reports counts per patch

### 5. Gaia Validation

Implemented in `gaia_validation.py`. This stage:

- Determines the relevant sky footprint and queries Gaia sources within it using `astroquery`
- Cross-matches SExtractor detections with Gaia sources and calculates source separations
- Produces a combined catalog with columns: `NUMBER, MAG_AUTO, MAGERR_AUTO, X_IMAGE, Y_IMAGE, FLAGS, ELLIPTICITY, CLASS_STAR, FLUX_RADIUS, RA, Dec, MAG_CAL, GAIA_RA, GAIA_Dec, GAIA_Gmag, Sep_arcsec, TILE, SOURCE, REGION`

The Gaia G-band magnitude is used in the subsequent stellar-selection stage.

### 6. Stellar Selection and Patch Definition

Implemented in `star_selection.py`. This stage:

- Applies selection criteria of `CLASS_STAR > 0.8` and `classprob_dsc_combmod_star > 0.98` to identify high-confidence stars
- Overlays confirmed stars on the imaging data for visual validation
- Generates the magnitude-vs-flux-radius plot (flux radius from the SExtractor catalog, Gaia G magnitude as stellar magnitude)
- Divides selected stars into spatial/morphological patches — core, intermediate, and optionally two outer regions — each associated with a particular cutout scale, controlled through `config.py`

This grouping allows the pipeline to use different stellar cutout sizes for different portions of the final radial PSF profile.

### 7. Cutout Generation

Implemented in `cutouts.py`. This stage:

- Redefines the four patch regions, each tagged with a cutout size (200–2000 px) and folder name, and creates one `cutouts/` output folder per patch
- `get_patch()`: given a star's flux radius and magnitude, returns which patch (and cutout size) it belongs to
- `is_valid_cutout()`: rejects cutouts where more than 20% of pixels are NaN/zero
- `generate_cutout()`: for each tile, loads the tile image and its matched Gaia catalog; for each star, determines its patch, checks it isn't too close to the image edge, extracts a `Cutout2D`, validates it, and saves it as a FITS file named by tile/star ID/size; tracks and prints save/skip counts per tile and totals
- For each patch folder, randomly samples 5 saved cutouts and displays them in a row with percentile scaling for quality control

### 8. Contamination Masking

Implemented in `masking.py`, which imports `cutouts.generate_cutout()` rather than re-extracting cutouts. Neighboring astronomical sources can contaminate the extended wings of a stellar PSF; PoLaRiS uses SExtractor segmentation maps to identify these. For each stellar cutout, the stage:

- Loads the cutout and the matching segmap cutout
- Identifies the target object ID at the cutout's center pixel; skips (and deletes files) if the center is heavily displaced
- Builds a primary mask (star+background pixels = 0, rest = 1), flags segmented pixels belonging to other objects, and dilates that mask by 4 iterations using a full 8-connectivity structuring element to grow contamination flags outward
- Builds a final mask as the inverse of the primary mask, and computes the multiplied image (cutout × final mask) to zero out contaminating neighbors
- Saves all five products per star: cutout, segmap cutout, primary mask, final mask, and multiplied image, organized into `cutouts/`, `segmaps/`, `primary_masking/`, `final_masking/`, and `multiplied/`
- Preview: for 5 random stars per patch, shows a 5-panel figure (Cutout, Segmap, Primary Mask, Final Mask, Masked Image)

### 9. PSF Stacking

Implemented in `stacking.py`. This stage:

- Loads masked stellar cutouts via `load_masked_cutout()`, representing invalid/masked pixels as NaN
- Applies quality-rejection criteria — maximum masked fraction, expected cutout dimensions, and centering tolerance (via `peak_offset()`) — configured per patch in `config.py`
- Normalizes each accepted cutout using the median flux within a configured annulus (`annulus_median()`) to establish a common flux scale before stacking
- Combines normalized cutouts to produce a median PSF stack, a mean PSF stack, and a per-pixel uncertainty map
- Reports acceptance/rejection statistics per patch and visualizes the stacks using log stretching and the normalization annulus (`plot_stack()`)

The median stack provides a robust representative PSF while reducing the influence of individual outliers.

### 10. Radial Profile Construction

Implemented in `radial_profile.py`. For each patch, the pipeline:

- Loads the median-stacked PSF FITS file
- Computes a radial profile using `photutils.RadialProfile` — a linear progression from 1–10 px and log-spaced radii from 10–1500 px — normalized to peak = 1
- Applies a configured vertical shift (`YSHIFT`) per patch to visually align overlapping profiles, and plots all patches on a hybrid plot
- Stitches the individual patch profiles into one continuous curve using configured boundaries — `CORE_MAX=40`, `INT_MAX=90`, `O1_MAX=245` — so that each patch contributes only over its designated radial interval

The radial sampling combines linear and logarithmic intervals, sampling the central PSF more densely while extending the profile to larger radii. The resulting profile is a single continuous representation of the PSF over the combined radial range covered by the individual stacks. The exact boundaries are configuration parameters and should be reviewed when applying PoLaRiS to a new dataset.

## Helper Functions

Common functionality shared across pipeline stages is implemented in `helpers.py`:

- **`make_radius_map()`** — builds a radial-distance-from-center grid
- **`annulus_median()`** — median flux within an annulus, used for cutout normalization
- **`peak_offset()`** — distance of the brightest pixel from the geometric center, used as a centering-quality criterion
- **`load_masked_cutout()`** — loads a cutout and applies its multiplied mask, converting masked/zero pixels to NaN
- **`get_patch()`** — determines the patch associated with a source based on the configured selection regions
- **`is_valid_cutout()`** — checks whether a stellar cutout satisfies the configured validity criteria
- **`plot_stack()`** — log-stretches and displays median/mean stacks side by side with colorbars and drawn normalization-annulus rings

## Pipeline Execution

The complete workflow is orchestrated by `run_pipeline.py`, which runs the ten stages above in order:

```bash
python run_pipeline.py
```

Before running the pipeline, review `polaris/config.py` and set the appropriate input/output paths, tile names, patch definitions, cutout sizes, normalization annuli, masking criteria, centering criteria, and profile-stitching boundaries.

## Example Notebooks

Example Jupyter notebooks are provided in `notebooks/` and walk through the pipeline interactively, stage by stage, using a small example dataset rather than a full survey tile. They are intended for exploration and for reviewers who want to inspect intermediate outputs (segmentation maps, cutout previews, mask panels, stacks, radial profiles) without running the full `run_pipeline.py` script end to end.

| Notebook | Covers |
|---|---|
| `01_preprocessing_and_detection.ipynb` | FITS inspection, downsampling, SExtractor source detection |
| `02_background_and_point_sources.ipynb` | Background subtraction, point-source identification |
| `03_gaia_validation_and_selection.ipynb` | Gaia cross-matching, stellar selection, patch definition |
| `04_cutouts_and_masking.ipynb` | Cutout generation, contamination masking |
| `05_stacking_and_radial_profile.ipynb` | PSF stacking, radial profile construction and stitching |

Run them from the repository root after installing the package (`pip install -e .`) so the `polaris` modules import correctly:

```bash
jupyter lab notebooks/
```

*(Adjust the notebook filenames/table above once the actual notebooks are added — this is a placeholder outline you can edit to match what you build.)*

## Output Products

Depending on the configured workflow, PoLaRiS produces:

- SExtractor catalogs and segmentation maps
- Background maps and background-subtracted images
- Gaia-matched and selected stellar-source catalogs
- Stellar FITS cutouts and segmentation-map cutouts
- Primary masks, final contamination masks, and masked/multiplied cutouts
- Median and mean PSF stacks, and per-pixel uncertainty maps
- Patch-specific radial profiles and the continuous stitched radial PSF profile

The exact directory organization and filenames are controlled by the pipeline configuration.

## Reproducibility and Configuration

PoLaRiS centralizes processing parameters in `config.py` rather than distributing them across individual modules. For reproducible processing, users should preserve, alongside their results:

- The PoLaRiS version or commit
- The input dataset and survey release
- The relevant configuration parameters
- The Python environment and dependency versions

When applying PoLaRiS to data from a different survey, the configuration should be reviewed rather than assuming that parameters used for KiDS or the DESI Legacy Imaging Surveys are directly transferable.

## Testing

Tests are provided in `tests/` and can be run with:

```bash
pytest
```

For a full scientific run, users should additionally inspect the diagnostic outputs generated at each stage — source detections, segmentation maps, Gaia cross-matches, stellar-selection plots, cutouts, contamination masks, stacked PSFs, and radial profiles — since these visual checks aren't something an automated test can fully substitute for.

## Current Scope and Limitations

- The current configuration is optimized for KiDS and DESI Legacy Imaging Survey imaging
- SExtractor is required for the current source-detection and segmentation workflow
- Gaia cross-matching is currently part of the stellar-validation workflow
- Patch definitions, cutout sizes, and profile-stitching boundaries are configuration-dependent

Applying the pipeline to another survey may require changes to image preprocessing, source-selection criteria, photometric quantities, pixel scale, patch definitions, and other configuration parameters. The software should not be assumed to be survey-independent without appropriate validation.

## Reference

| # | Title | Link |
|---|---|---|
| 1 | The Buildup of the Intracluster Light of A85 as Seen by Subaru's Hyper Suprime-Cam | IOP Science |
| 2 | The Hyper Suprime-Cam Extended Point Spread Functions and Applications | MNRAS |

Additional methodological references should be cited in the associated JOSS paper where specific algorithms or scientific procedures are discussed.

### Software and Data Resources

PoLaRiS makes use of the following astronomical software and data resources: SExtractor for source detection and segmentation-map generation; Astropy for astronomical data structures, FITS handling, coordinates, and image operations; Photutils for image and photometric analysis, including radial-profile construction; Astroquery for querying astronomical archives; Gaia for stellar-source validation; and imaging from the Kilo-Degree Survey (KiDS) and the DESI Legacy Imaging Surveys.

## Acknowledgements

The authors acknowledge the developers and maintainers of the open-source astronomical software used by PoLaRiS, including SExtractor, Astropy, Photutils, and Astroquery, and the teams responsible for the Gaia, KiDS, and DESI Legacy Imaging Surveys data products used during development and testing.

- SExtractor: Bertin, E. & Arnouts, S. (1996), *Astronomy and Astrophysics Supplement Series*, 117, 393–404
- KiDS Survey: https://kids.strw.leidenuniv.nl/
- DESI Legacy Surveys: https://www.legacysurvey.org/

Last Updated: 19 August, 2026
---

*Pipeline documentation for PSF Modeling.*
