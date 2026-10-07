# Marine SDM Benchmarking

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23202487.svg)](https://doi.org/10.5281/zenodo.23202487)

Benchmark and harmonisation protocol for marine species distribution models under data limitation. Seven methods, four FISHGLOB surveys, factorial design over sample size, spatial coverage and inference target; validation optimism, metric translation and the `smdMetricHarmonise` R package. One command, resumable, multi core execution in R.

This repository accompanies the manuscript *A benchmark harmonisation protocol for species distribution modelling under data limitation* (Afrifa-Yamoah, in preparation). It holds the complete pipeline, the full scale results and the pilot results, the theoretical verification, the tables and figures of the Supporting Information, and the package; the manuscript itself is not distributed here.

Author: Ebenezer Afrifa-Yamoah, School of Science, Edith Cowan University, ORCID [0000-0003-1741-9249](https://orcid.org/0000-0003-1741-9249).

## What the study does

Species distribution models are compared on performance figures produced under different sample sizes, spatial designs, validation schemes and metrics, which makes those figures hard to reconcile. The benchmark holds the data generating process fixed, by drawing all training and test sets from the same real bottom trawl surveys, and varies the confounds factorially:

| Factor | Levels |
|---|---|
| Surveys | Eastern Bering Sea, North Sea (IBTS), Gulf of Mexico, Northeast US shelf (FISHGLOB) |
| Species per survey | 15 at full scale (5 rare, 5 intermediate, 5 common); 3 in the pilot |
| Training sample size | 30, 50, 100, 200, 500 |
| Spatial coverage | clustered, intermediate, distributed |
| Inference target | conditional (interpolation), marginal (spatial extrapolation to a held-out block) |
| Methods | RF, spatial RF, BRT, presence only MaxEnt, spatial GAM, geostatistical GP, joint SDM |
| Replicates | 100 at full scale; 4 in the pilot |

Each fit is scored on nine metrics. The analysis quantifies how performance scales with sample size, the effect of spatial structure relative to algorithm choice, the coverage trade off, validation optimism and its decomposition, and learns a cross-metric harmonisation protocol with an explicit reliability class for every translation. A theoretical foundation (manuscript Section 3 and Appendices A to C) proves which translations are identities, why discrimination cannot be translated to error, the cost of pooling heterogeneous surveys, and why the marginal performance of a spatial predictor is not identifiable from interpolation data; `R/10_theory_analysis.R` checks those results against the benchmark's own scores.

## Repository layout

```
run_full_local.R            one command full run (see docs/full_scale_execution.md)
R/
  config.R                  design, scale switch (pilot | full), paths, core count
  00_download_data.R        fetch the four FISHGLOB surveys into data_raw/
  01_preprocess.R           modelling frames and species selection
  02_subsample.R            factorial subsampler, held-out blocks, draw geometry
  03_methods.R              the seven method wrappers
  04_metrics.R              performance metrics
  05_helpers.R              shared cell logic (same fit contrast, isotonic REF and CAL)
  05_run_sim.R              serial driver (pilot)
  05_run_sim_parallel.R     resumable, self healing parallel driver (full scale)
  05_run_sim_shard.R        socket free serial driver for one shard of the design (clusters, sandboxes)
  05c_merge_parts.R         merge replicate block checkpoints into results/fit_metrics.rds
  05b_oracle.R              near oracle intrinsic term and variogram ranges
  06_analysis.R             hierarchical models, variance partition, optimism
  06b_theory_full.R         analyses that need the full scale quantities
  07_harmonise.R            metric translations, reliability, package assembly
  08_figures.R              result figures
  09_explore_spatial.R      sampling profile, occurrence and design figures
  07b_package_theory.R      refresh the package's theory table once the bootstrap of 10 is complete
  10_theory_analysis.R      identity verification, bootstrap, theory tests
  11_manuscript_numbers.R   every number cited in the manuscript, with its definition and source
  12_supplement.R           tables S1 to S16 and figures S1 to S9 of the Supporting Information
  palette.R                 shared colour vocabulary: a hue identifies a method and nothing else
data_processed/             full scale modelling frames (regenerable from 01)
results/full/               full scale tables, manuscript_numbers.csv, bootstrap classes, near oracle fits
results/full/supplement/    Supporting Information tables (table_S*.csv)
figures/full/               full scale manuscript figures; figures/full/supplement/ holds figures S1 to S9
results/pilot/              every table reported in the pilot manuscript
figures/pilot/              every figure reported in the pilot manuscript
smdMetricHarmonise/         R package source (version 0.3.0)
docs/full_scale_execution.md   requirements, runtime, resuming, troubleshooting
```

Raw survey files are not committed; `R/00_download_data.R` fetches them (41 MB) from the public FISHGLOB repository, and `run_full_local.R` does so automatically when they are absent.

## Quick start

```r
install.packages(c("data.table", "mgcv", "ranger", "gbm", "glmnet", "fields",
                   "lme4", "pROC", "ggplot2", "viridis", "patchwork", "maps",
                   "MASS", "boot"))
```

Reproduce the pilot (about 10 minutes on one core):

```bash
SMD_SCALE=pilot Rscript R/01_preprocess.R
SMD_SCALE=pilot Rscript R/05_run_sim.R
Rscript R/06_analysis.R; Rscript R/07_harmonise.R; Rscript R/08_figures.R
```

Run the full study (1,260,000 fits, roughly 60 to 120 core hours):

```bash
SMD_CORES=12 Rscript run_full_local.R
```

The run checkpoints every design cell, so it can be stopped and resumed, split across machines by survey, and reduced for a first pass with `SMD_REPLICATES=25`. Details, runtime estimates and troubleshooting are in `docs/full_scale_execution.md`.

## Full scale results

The full study planned 1,260,000 fits (4 surveys × 15 species × 30 design cells × 100 replicates × 7 methods). Of 180,000 training and test splits, 11,747 were degenerate (all among rare and intermediate species), leaving 1,177,771 scored fits, of which 1,177,622 converged with finite metrics. Every table and figure of the manuscript and its Supporting Information is in `results/full/` and `figures/full/`, and `results/full/manuscript_numbers.csv` lists every number cited in the text with its definition and the file it comes from.

After `run_full_local.R` has finished, the remaining stages are:

```bash
export SMD_SCALE=full SMD_REPLICATES=100
Rscript R/10_theory_analysis.R      # repeat until it reports 100 of 100 bootstrap replicates
Rscript R/07b_package_theory.R      # bring the package's theory table up to date
Rscript R/11_manuscript_numbers.R   # results/manuscript_numbers.csv
Rscript R/12_supplement.R           # results/supplement/ and figures/supplement/
```

The fit level table `results/fit_metrics.rds` (1,177,771 rows, about 150 MB) is not committed because of its size; it will be archived with a DOI, and every script above reads it from `results/`.

## The harmonisation package

```r
install.packages("smdMetricHarmonise", repos = NULL, type = "source")  # from a built tarball
library(smdMetricHarmonise)
harmonise_metric(0.75, from = "AUC", to = "TSS")
harmonise_metric(0.10, from = "R2", to = "RMSE", prevalence = 0.25)   # exact by identity
```

Translations among R2, RMSE and the Brier score are returned by identity when the evaluation prevalence is supplied. Learned translations return the estimate, an approximate one standard deviation band, their reliability class, the stability of that class under a cluster bootstrap over species, and the binding source of instability (inference target or survey); a request for an unreliable translation raises a warning. Build with `R CMD build smdMetricHarmonise` after `R/07_harmonise.R` has run.

## Reproducibility

Release v1.0.0 of this repository is archived on Zenodo, https://doi.org/10.5281/zenodo.23202487, and is the version used for the results reported in the manuscript. A global seed drives a deterministic seed for every design cell, so each fit is reproducible regardless of the number of cores or the order in which cells run. Checkpoint folders are keyed to scale and replicate count, so pilot and full scale results never mix. The full scale tables and figures were produced from 1,177,771 scored fits with R 4.3.3. The pilot tables and figures in this repository were produced from 9,800 fits (4 surveys × 3 species × 30 design cells × 4 replicates × 7 methods, less 40 degenerate splits).

## Data

Survey data are from FISHGLOB (Maureaud et al. 2024, *Scientific Data*), https://github.com/fishglob/FishGlob_data, used under its licence.

## How to cite

Afrifa-Yamoah, E. (2026). *Marine SDM Benchmarking: benchmark harmonisation protocol for species distribution models under data limitation* (Version 1.0.0) [Computer software]. Zenodo. https://doi.org/10.5281/zenodo.23202487

Please also cite the accompanying article once it is published. GitHub's "Cite this repository" button gives the same reference from `CITATION.cff`.

## Licence

MIT. See `LICENSE`.
