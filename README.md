# PoLaRiS: Point-spread Function Reconstruction of Stellar Sources

**PoLaRiS** is an automated pipeline for point spread function (PSF) modeling and reconstruction from multi-band astronomical imaging data. It is currently optimized for the **Kilo-Degree Survey (KiDS)** and the **DESI Legacy Imaging Surveys**.

---

## Table of Contents

- [Step 0: Getting Started](#step-0-getting-started)
- [Step 1: Configuration](#step-1-configuration)
- [Step 2: Data Preprocessing](#step-2-data-preprocessing)
- [Step 3: Background Subtraction](#step-3-background-subtraction)
- [Step 4: Identify Point Sources](#step-4-identify-point-sources)
- [Step 5: MAG_AUTO vs FLUX_RADIUS Plot](#step-5-mag_auto-vs-flux_radius-plot)
- [Step 6: Gaia Validation](#step-6-gaia-validation)
- [Step 7: Star Selection](#step-7-star-selection)
- [Step 8: Cutout Generation](#step-8-cutout-generation)
- [Step 9: Masking](#step-9-masking)
- [Step 10: Stacking](#step-10-stacking)
- [Step 11: Radial Profile Construction](#step-11-radial-profile-construction)
- [utils.py](#utilspy)
- [References](#references)
- [Acknowledgements](#acknowledgements)

---

## Pipeline Overview

```
Raw tile → SExtractor → Background subtraction → Point-source selection
   → Gaia validation → Star grouping (Core / Intermediate / Outer_1 / Outer_2)
   → Cutouts → Masking → Normalized median stacks → Stitched radial profile
```

---

## Step 0: Getting Started

This step sets up all required dependencies for PSF modeling, including both the Python libraries and SExtractor.

### 0.1 Clone the repository

```bash
git clone https://github.com/Radit-Raian/PoLaRiS.git
cd PoLaRiS
pip install -r requirements.txt
```

### 0.2 Build SExtractor

```bash
git clone https://github.com/astromatic/sextractor
cd sextractor
sh autogen.sh
./configure
make -j
sudo make install
sex   # verify the installation
```

---

## Step 1: Configuration

Modify the SExtractor configuration template as needed. Additional parameters are controlled later through this configuration file.

---

## Step 2: Data Preprocessing

In this step, we will:

- Visualize the image and FITS header to verify data quality
- Downsample the FITS image for a quick quality inspection
- Set the pixel scale from the survey
- Define input/output directories for the SExtractor work
- Point to the necessary configuration files (`.nnw`, `.psf`, `.param`, convolution filter, etc.)
- Write a per-tile SExtractor configuration file and loop through all tiles
- Run SExtractor using the current configuration template to detect primary sources
- Visualize the downsampled segmentation map to validate source detection (KiDS FITS files are large, so downsampling avoids desktop crashes)

---

## Step 3: Background Subtraction

In this step, we will:

- For each tile, load the raw image, the SExtractor background map, and the background RMS map
- Compute `background-subtracted = raw − background`
- Generate a background-subtracted image for further analysis
- Perform background subtraction using Photutils

---

## Step 4: Identify Point Sources

In this step, we will:

- For the same tiles, parse the SExtractor `.cat` header lines to reconstruct the column names
- Load the catalog into a DataFrame and filter for point sources (`CLASS_STAR ≥ 0.84`)
- Load the downsampled tile image and overlay red circle markers at the (scaled) pixel positions of the identified point sources

---

## Step 5: MAG_AUTO vs FLUX_RADIUS Plot

In this step, we will:

- Compute an adaptive `CLASS_STAR` threshold per tile by inspecting the plot
- Distribute the selected stars into four groups based on their position on the plot: **Core**, **Intermediate**, **Outer_1**, and **Outer_2**
- Clip the selected stars to acceptable flux-radius and magnitude ranges for cleaner analysis
- Generate counts per patch to get the total number of stars in each group
- Optionally define two patches under the same group to refine the selection

---

## Step 6: Gaia Validation

In this step, we will:

- Query Gaia sources within the image footprint using `astroquery` over the defined rectangular sky region
- Cross-match SExtractor-detected sources with the Gaia catalog
- Validate the star magnitudes against the obtained Gaia data
- Write a `.cat` file with the following columns:

```
NUMBER, MAG_AUTO, MAGERR_AUTO, X_IMAGE, Y_IMAGE, FLAGS, ELLIPTICITY, CLASS_STAR,
FLUX_RADIUS, RA, Dec, MAG_CAL, GAIA_RA, GAIA_Dec, GAIA_Gmag, Sep_arcsec,
TILE, SOURCE, REGION
```

---

## Step 7: Star Selection

In this step, we will:

- Apply the selection criteria `CLASS_STAR > 0.8` and `classprob_dsc_combmod_star > 0.98` to identify high-confidence stars
- Overlay confirmed stars on the imaging data for visual validation
- Generate the **Magnitude vs. Flux Radius 2.0** plot, with the SExtractor `FLUX_RADIUS` on the x-axis and the Gaia G-band magnitude on the y-axis
- Define patches that divide the PSF into core, intermediate, and (optionally) two outer regions

---

## Step 8: Cutout Generation

In this step, we will:

- Redefine the same four patch regions, now each tagged with a cutout size (200–2000 px) and a folder name
- Create one `cutouts/` output folder per unique patch name
- Use the helper `get_patch()`: given a star's flux radius and magnitude, it returns the patch (and cutout size) the star belongs to
- Use the helper `is_valid_cutout()`: rejects cutouts where more than 20% of pixels are NaN/zero
- Loop over all tiles, loading each tile image and its matched Gaia catalog. For each star we:
  - determine its patch,
  - check that it is not too close to the image edge,
  - extract a `Cutout2D`,
  - validate it, and
  - save it as a FITS file named by tile / star ID / size
- Track and print save/skip counts per tile and in total
- For each patch folder, randomly sample 5 saved cutouts and display them in a row with percentile scaling

---

## Step 9: Masking

In this step, we will:

- Add new subfolders per patch: `cutouts`, `segmaps`, `primary_masking`, `final_masking`, `multiplied`
- For each star, extract a cutout and the matching segmentation-map cutout
- Identify the "target" object ID at the cutout's center pixel; skip (and delete files) if the center is heavily displaced
- Build a **primary mask** by writing the star + background pixels as `0` and the rest as `1`
- Flag all segmented pixels belonging to other objects (not the central star), then dilate that mask by 4 iterations using a full 8-connectivity structuring element to grow contamination flags outward
- Build a **final mask** as the inverse of the primary mask
- Compute the **multiplied image** (`cutout × final mask`) to zero out contaminating neighbors
- Save all five products per star: cutout, segmap cutout, primary mask, final mask, multiplied image
- Preview: for 5 random stars per patch, show a 5-panel figure (Cutout, Segmap, Primary Mask, Final Mask, Masked Image)

---

## Step 10: Stacking

In this step, we will define per-patch parameters: normalization annulus radii, maximum allowed masked fraction, expected cutout size, centering tolerance, and stack grouping (currently each patch stacks alone).

**Helper functions**

| Function | Purpose |
|---|---|
| `make_radius_map()` | Builds a radial-distance-from-center grid |
| `annulus_median()` | Median flux within an annulus (used for normalization) |
| `peak_offset()` | Distance of the brightest pixel from the geometric center (used to reject miscentered stars) |
| `load_masked_cutout()` | Loads a cutout and applies its multiplied mask, converting masked/zero pixels to NaN |
| `plot_stack()` | Log-stretches and displays median/mean stacks side by side with colorbars and the drawn normalization-annulus rings |

**Procedure**

- Load every cutout and apply the rejection criteria in sequence
- Normalize accepted cutouts by their annulus median flux and collect them into a stack
- Compute the median stack, the mean stack, and a per-pixel uncertainty map
- Plot the median/mean stacks
- Print a summary table of acceptance/rejection counts per patch

---

## Step 11: Radial Profile Construction

In this step, we will:

- Load each patch's median-stacked PSF FITS file, e.g. `Core/stacked/median_psf_Core.fits`
- Compute a radial profile using `photutils.profiles.RadialProfile` (linear spacing from 1–10 px and log-spaced radii from 10–1500 px), normalized to peak = 1
- Apply a manual y-axis shift per patch (to visually align the overlapping profiles) and plot all four patches on a hybrid plot
- Define stitching boundaries to combine the four separate radial profiles into one continuous curve, using each patch only within its designated radius range:

| Patch | Radius range (px) |
|---|---|
| Core | ≤ `CORE_MAX` = 40 |
| Intermediate | 40 – `INT_MAX` = 90 |
| Outer_1 | 90 – `O1_MAX` = 250 |
| Outer_2 | > 250 |

- Plot the final continuous PSF radial profile, with the bottom axis in pixels and a secondary top axis in degrees (using `pixscale = 0.2` arcsec/pixel, i.e. `deg = px × pixscale / 3600`)

> **Note:** adjust `CORE_MAX`, `INT_MAX`, `O1_MAX` and the per-patch y-shifts to match where the profiles agree for your data.

---

## utils.py

Shared helper functions used across the pipeline steps.

---

## References

| # | Title | Link |
|---|---|---|
| 1 | The Buildup of the Intracluster Light of A85 as Seen by Subaru's Hyper Suprime-Cam | IOP Science |
| 2 | The Hyper Suprime-Cam extended point spread functions and applications | MNRAS |

---

## Acknowledgements

- **SExtractor:** Bertin, E. & Arnouts, S. (1996), *Astronomy and Astrophysics Supplement Series*, 117, 393–404
- **KiDS Survey:** <https://kids.strw.leidenuniv.nl/>
- **DESI Legacy Surveys:** <https://www.legacysurvey.org/>

---

*Pipeline documentation for PSF modeling. Last updated: 2026-04-29*
