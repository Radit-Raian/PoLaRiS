# PoLaRiS: Point-spread Light Reconstructions for Stars

PoLaRiS (Point-spread Light Reconstruction for Stars) is an automated pipeline for point spread function (PSF) modeling and reconstruction from multi-band astronomical imaging data. It is currently optimized for the Kilo-Degree Survey and the DESI Legacy Imaging Surveys.

## Table of Contents

1. [Getting Started](#step-0-getting-started)
2. [Configuration](#step-1-configuration)
3. [Data Preprocessing](#step-2-data-preprocessing)
4. [Background Substraction](#step-3-background-substraction)
5. [GAIA Validation](#step-4-gaia-validation)
6. [Star Selection](#step-5-star-selection)
7. [Cutout Generate](#step-6-cutout-generate)
8. [Masking](#step-7-masking)
9. [Stacking](#step-8-stacking)
10. [Radial Profile](#step-9-radial-profile)
11. [Reference](#references)
12. [Acknowledgements](#acknowledgements)

---

## Step 0: Getting Started

This step sets up all required dependencies for PSF modeling, including both Python libraries and SExtractor. 

### 0.1 Clone the repository

```bash
git clone https://github.com/Radit-Raian/PoLaRiS.git
cd PoLaRiS
pip install -r requirements.txt
```
### 0.2 Building SExtractor

```bash
git clone https://github.com/astromatic/sextractor
cd sextractor
sh autogen.sh
./configure
make -j
sudo make install
sex
```

## Step 1: Configuration

## Step 2: Data Preprocessing

## Step 3: Background Substraction

## Step 4: GAIA Validation

## Step 5: Star selection

## Step 6: Cutout Generate

## Step 7: Masking

## Step 8: Stacking

## Step 9: Radial Profile




## References

| # | Title | Link |
|---|-------|------|
| 1 | The Buildup of the Intracluster Light of A85 as Seen by Subaru's Hyper Suprime-Cam | [IOP Science](https://iopscience.iop.org/article/10.3847/1538-4357/abddb6) |
| 2 | The Hyper Suprime-Cam extended point spread functions and applications | [MNRAS](https://academic.oup.com/mnras/article/531/2/2517/7676869) |


---

## Acknowledgements

- **SExtractor:** Bertin, E. & Arnouts, S. (1996), Astronomy and Astrophysics Supplement Series, 117, 393-404
- **KiDS Survey:** https://kids.strw.leidenuniv.nl/ 
- **DESI Legacy Surveys:** https://www.legacysurvey.org/


---
  
*Pipeline documentation for PSF Modeling*  
*Last updated: 2026-29-04*
