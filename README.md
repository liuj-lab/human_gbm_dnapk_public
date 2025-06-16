# human_gbm_dnapk
A repository for analysis relating to human GBM models and DNA-PK

## Pre-commit
This repository uses pre-commit hooks to style and maintain standards for Python and R code. This requires
contributors to have [pre-commit](https://pre-commit.com/) installed. Upon installation, run `pre-commit install`
to automatically install the pre-commit hooks specified in `.pre-commit-config.yaml`.

Note that if running on C4, you'll need to have a version of `RScript` in the system's `PATH`. Running `module load CBI rstudio-server-controller`
does not do this automatically. You'll need to run `module load CBI r`.
