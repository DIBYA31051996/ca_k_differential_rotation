# KoSO Ca II K Differential Rotation

Reproducible analysis and figure-generation code for **“Differential Rotation of the Solar Chromosphere: A Century-long Perspective from Kodaikanal Solar Observatory Ca II K Data.”**

**Authors:** Dibya Kirti Mishra, Srinjana Routh, Bibhuti Kumar Jha, Theodosios Chatzistergos, Judhajeet Basu, Subhamoy Chatterjee, Dipankar Banerjee, and Ilaria Ermolli

[![DOI](https://img.shields.io/badge/DOI-10.3847%2F1538--4357%2Fad1188-blue)](https://doi.org/10.3847/1538-4357/ad1188)

## Overview

This repository contains the IDL analysis routines, processed data products, and Python notebooks used to study chromospheric differential rotation from the calibrated Kodaikanal Solar Observatory (KoSO) Ca II K archive spanning 1907–2007.

The analysis measures the displacement of chromospheric structures between observations through image cross-correlation. The resulting mean rotation profile is

$$
\Omega(\theta) = 14.61 \pm 0.04 - (2.18 \pm 0.37)\sin^2\theta - (1.10 \pm 0.61)\sin^4\theta
\quad [^\circ\,\mathrm{day}^{-1}].
$$

The inferred equatorial chromospheric rotation is 1.59% faster than the photospheric rate derived from sunspots. The study finds no significant north–south asymmetry or systematic century-scale variation.

> [!IMPORTANT]
> **Main KoSO image-correlation code:** [`src/idl/koso_image_correlation.pro`](src/idl/koso_image_correlation.pro)
>
> This is the primary implementation of the measurement workflow. It reads consecutive observations, corrects image orientation, masks the off-disk region, remaps the images into heliographic coordinates, divides them into 5° latitude strips, cross-correlates matching strips, rejects unreliable peaks, and converts the measured longitudinal shifts into sidereal angular rotation rates.

## Image-correlation workflow

1. Read and orient two consecutive full-disk observations.
2. Mask pixels outside the adopted solar radius.
3. Transform each disk from Cartesian image coordinates to heliographic longitude and latitude.
4. Extract matching 5° latitude bands over ±55° longitude.
5. Estimate the expected longitudinal offset from a reference differential-rotation law.
6. Cross-correlate each pair of strips and locate the correlation maximum.
7. Reject weak or poorly defined correlation peaks.
8. Divide the accepted angular displacement by the observation interval and apply the sidereal correction.
9. Average accepted measurements by latitude and fit the differential-rotation profile.

Supporting IDL routines—including coordinate conversion, disk masking, time conversion, sidereal correction, and profile fitting—are located in [`src/idl/`](src/idl/). The separate [`cross_correlation2.pro`](src/idl/cross_correlation2.pro) program applies the correlation approach to MDI data for method validation.

## Repository layout

```text
.
├── config/                 Matplotlib style configuration
├── data/                   Processed IDL save files used by the notebooks
│   └── legacy/             Legacy text outputs and auxiliary analysis data
├── docs/images/            Images embedded in the documentation
├── figures/                Generated and published PDF figures
├── notebooks/              Python figure-reproduction notebooks
├── src/idl/                Main IDL analysis code and supporting routines
├── CITATION.cff            Citation metadata
├── LICENSE                 Project license
└── requirements.txt        Python dependencies
```

## Figure notebooks

| Notebook | Description |
| :--- | :--- |
| [`notebooks/method_figures.ipynb`](notebooks/method_figures.ipynb) | Reproduces the methodology, image-remapping, cross-correlation, and quality-threshold figures. |
| [`notebooks/result_figures.ipynb`](notebooks/result_figures.ipynb) | Reproduces the differential-rotation profiles, validation comparisons, and temporal results. |

The notebooks read processed `.sav` files from `data/`, use the shared style in `config/`, and write PDFs to `figures/`.

## Average chromospheric rotation profile

![Average chromospheric differential-rotation profiles from KoSO Ca II K observations and comparisons with other measurements and observatories](docs/images/figure4_main_result.png)

Figure 4 of the paper. Panel (a) shows the average KoSO Ca II K rotation profile for 1907–2007 and the profiles before and after 1980. Panel (b) compares the result with selected photospheric and chromospheric measurements from the literature. Panel (c) compares KoSO, Rome/PSPT, and Meudon Ca II K measurements over 2000–2002.

## Running the figure notebooks

Create an isolated Python environment and install the dependencies:

```bash
python -m venv .venv
```

On Linux or macOS:

```bash
source .venv/bin/activate
python -m pip install -r requirements.txt
jupyter lab notebooks/
```

On Windows PowerShell:

```powershell
.venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
jupyter lab notebooks/
```

Run notebook cells from the `notebooks/` directory so their `../data/`, `../config/`, and `../figures/` paths resolve correctly.

The legacy analysis is written in IDL and relies on astronomy/SolarSoft routines such as `READ_SDO`, `CORREL_IMAGES`, and `CORRMAT_ANALYZE`. Adjust the input directory at the beginning of the main IDL program before running it in a configured IDL environment.

The analysis files are callable IDL routines. For example:

```idl
KOSO_IMAGE_CORRELATION
DIFFERENTIA_ASYMMETRY
CYCLESTRENGTH
```

IDL automatically compiles a routine when its `.pro` filename matches the routine name and `src/idl/` is included in `!PATH`.

## Data availability

The original analysis workspace did not contain `south_hem.sav`, `wilson_scatter_new.sav`, or `wilson_scatter_new1.sav`. Cells that use these processed inputs cannot be regenerated until the files are restored; the corresponding published PDFs remain available in `figures/`. Original observatory data are available from the archives cited in the paper.

## Citation

If you use this repository, cite the associated paper using the metadata in [`CITATION.cff`](CITATION.cff):

> Mishra, D. K., et al. (2023). *Differential Rotation of the Solar Chromosphere: A Century-long Perspective from Kodaikanal Solar Observatory Ca II K Data*. The Astrophysical Journal. https://doi.org/10.3847/1538-4357/ad1188

## License

See [`LICENSE`](LICENSE) for licensing terms.
