# PoLaRiS Demo Notebooks

Demo notebook(s) showing how to run the [PoLaRiS](https://github.com/Radit-Raian/PoLaRiS)
pipeline end to end. Each cell calls one pipeline step's function directly —
no reimplementation, just the same functions the CLI (`python -m polaris.<step>`)
calls under the hood.

## Prerequisites

- `polaris` installed and importable (`pip install -e .` from the PoLaRiS repo,
  or add it to your `PYTHONPATH`)
- SExtractor installed and on `PATH`
- Your own FITS tile(s) already in `input/FITS_images/`, and the SExtractor
  support files already in `input/sextractor/`

## Contents

| Notebook | Covers |
|---|---|
| `polaris_demo.ipynb` | Preprocessing → SExtractor → background → point sources → Gaia validation → star selection → cutouts → masking → stacking → radial profile |

Before running, edit the `config.TILE_NAMES` line in the setup cell to point
at your tile(s).
