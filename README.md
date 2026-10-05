# Marine SDM Benchmarking

Benchmark and harmonisation protocol for marine species distribution models under data limitation. Seven methods, four FISHGLOB surveys, factorial design over sample size, spatial coverage and inference target; validation optimism, metric translation and the `smdMetricHarmonise` R package. One command, resumable, multi core execution in R.

This repository accompanies the manuscript *A benchmark harmonisation protocol for species distribution modelling under data limitation* (Afrifa-Yamoah, in preparation). It holds the complete pipeline, the pilot results reported in the manuscript, the theoretical verification, and the package; the manuscript itself is not distributed here.

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
  10_theory_analysis.R      identity verification, bootstrap, theory tests
data_processed/             full scale modelling frames (regenerable from 01)
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

## The harmonisation package

```r
install.packages("smdMetricHarmonise", repos = NULL, type = "source")  # from a built tarball
library(smdMetricHarmonise)
harmonise_metric(0.75, from = "AUC", to = "TSS")
harmonise_metric(0.10, from = "R2", to = "RMSE", prevalence = 0.25)   # exact by identity
```

Translations among R2, RMSE and the Brier score are returned by identity when the evaluation prevalence is supplied. Learned translations return the estimate, an approximate one standard deviation band, their reliability class, the stability of that class under a cluster bootstrap over species, and the binding source of instability (inference target or survey); a request for an unreliable translation raises a warning. Build with `R CMD build smdMetricHarmonise` after `R/07_harmonise.R` has run.

## Reproducibility

A global seed drives a deterministic seed for every design cell, so each fit is reproducible regardless of the number of cores or the order in which cells run. Checkpoint folders are keyed to scale and replicate count, so pilot and full scale results never mix. The pilot tables and figures in this repository were produced from 9,800 fits (4 surveys × 3 species × 30 design cells × 4 replicates × 7 methods, less 40 degenerate splits).

## Data

Survey data are from FISHGLOB (Maureaud et al. 2024, *Scientific Data*), https://github.com/fishglob/FishGlob_data, used under its licence.

## Licence

MIT. See `LICENSE`.
