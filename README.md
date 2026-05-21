# human_gbm_dnapk

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20334626.svg)](https://doi.org/10.5281/zenodo.20334626)

A repository for analysis relating to human GBM models and DNA-PKcs. Software dependencies and versions are indicated in the individual analysis scripts. These analyses can be run on normal desktop hardware with a typical install and run time of under 1 hour.

## Environment setup

Two environments are required to run code in this repository: an R environment managed with `renv`, and a Python environment managed with `conda`.

### R environment (renv)

R dependencies and their pinned versions are captured in `renv.lock`. To restore the exact environment used for analysis:

1. Install R 4.4.1 (the version pinned in `renv.lock`).
2. From the repository root, launch R and run:
   ```r
   install.packages("renv")   # if not already installed
   renv::restore()
   ```
   `renv::restore()` reads `renv.lock` and installs each package at the version it was originally run with. The activation script in `renv/activate.R` is sourced automatically through `.Rprofile`, so subsequent R sessions started from the repository root will use this project-local library.
3. Open `human_gbm_dnapk.Rproj` in RStudio (or run scripts via `Rscript` from the project root) to use the restored environment.

### Python environment (conda)

Python dependencies are captured in `environment.yml`. To create the environment:

```bash
conda env create -f environment.yml
conda activate human_gbm_dnapk
```

This creates a conda environment named `human_gbm_dnapk` with all required packages (scanpy, anndata, D-SPIN, MAGIC, scprep, etc.) at the versions used for the analyses. Launch Jupyter from within this environment to run any of the notebooks:

```bash
conda activate human_gbm_dnapk
jupyter lab
```

## Running an example

The full D-SPIN analysis on the GBM43 Perturb-seq dataset is reproducible end-to-end using two notebooks under `python/gbm43/`:

1. `dspin_preprocessing.ipynb` — QC, normalization, highly-variable-gene selection, PCA/UMAP/Leiden clustering. Produces the processed AnnData object consumed by the next step.
2. `dspin_run.ipynb` — oNMF gene-program discovery, BIC/MAGIC program-number selection, D-SPIN model fitting, GPT-4o-based program annotation, and downstream visualization.

### Download the example data

The raw input data and ancillary files required to run these two notebooks are archived on Zenodo:

> https://doi.org/10.5281/zenodo.20334626

Download the archive and place its contents under a top-level `data/` directory at the repository root. The notebooks expect, at minimum:

```
data/
├── gbm43_perturb_seq/
│   └── gbm43.h5ad
└── regev_lab_cell_cycle_genes.txt
```

`regev_lab_cell_cycle_genes.txt` is the standard S- and G2/M-phase marker list from Tirosh et al. 2016 (Regev lab). If it is not bundled with your Zenodo download, it can be fetched directly from the Theis lab's `scanpy_usage` repository:

```bash
curl -L -o data/regev_lab_cell_cycle_genes.txt \
  https://raw.githubusercontent.com/theislab/scanpy_usage/master/180209_cell_cycle/data/regev_lab_cell_cycle_genes.txt
```

(`dspin_preprocessing.ipynb` writes `gbm43_processed.h5ad` into `data/gbm43_perturb_seq/`, which is then loaded by `dspin_run.ipynb`.) Output figures and intermediate artifacts are written under a top-level `output/` directory, which will be created on first run.

### Run the notebooks

From the repository root, with the `human_gbm_dnapk` conda environment active:

```bash
cd python/gbm43
jupyter lab dspin_preprocessing.ipynb
# then, once preprocessing has produced gbm43_processed.h5ad:
jupyter lab dspin_run.ipynb
```

The GPT-4o annotation cells in `dspin_run.ipynb` require an OpenAI API key. Place it in a `.env` file at the repository root (e.g. `OPENAI_API_KEY=sk-...`); `python-dotenv` will pick it up.

### Other notebooks and scripts

All other notebooks under `python/` and scripts under `R/` and `shell/` are provided as reference and documentation for how the remaining analyses in the manuscript were run. Raw input objects required to re-execute them are not bundled in this repository, but can be provided upon request.

## Pre-commit
This repository uses pre-commit hooks to style and maintain standards for Python and R code. This requires
contributors to have [pre-commit](https://pre-commit.com/) installed. Upon installation, run `pre-commit install`
to automatically install the pre-commit hooks specified in `.pre-commit-config.yaml`.

Note that if running on C4, you'll need to have a version of `RScript` in the system's `PATH`. Running `module load CBI rstudio-server-controller`
does not do this automatically. You'll need to run `module load CBI r`.
