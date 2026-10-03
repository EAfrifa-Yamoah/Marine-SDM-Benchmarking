# Full scale local execution

The same code as the pilot, at the design scale: 4 surveys x 15 species x 30 design
cells x 7 methods x 100 replicates = **1,260,000 fits**, plus the near oracle stage.
Every stage is resumable; the run can be stopped and restarted at any time.

## Requirements
R >= 4.2 with: data.table, mgcv, ranger, gbm, glmnet, fields, lme4, pROC, ggplot2,
viridis, patchwork, maps, MASS, boot (all CRAN). `parallel` ships with R.

## Run
```bash
Rscript run_full_local.R                       # everything, all but one core
SMD_CORES=12 Rscript run_full_local.R          # choose cores
SMD_CORES=12 SMD_REPLICATES=25 Rscript run_full_local.R   # faster first pass
SMD_DATASETS=EBS Rscript run_full_local.R      # one survey at a time
```
Windows: set the variables with `$env:SMD_CORES=12` in PowerShell.

## Runtime
About 0.15 to 0.35 s per fit (the same fit contrast and isotonic recalibration add
roughly a third over the pilot), so 60 to 120 core hours for 100 replicates:
4 cores 15 to 30 h, 8 cores 8 to 15 h, 16 cores 4 to 8 h. A 25 replicate pass is one
quarter. The oracle stage is about 1,700 fits on up to 3,000 hauls each, roughly an
hour on one core; cap with SMD_ORACLE_MAX_N. Progress and ETA: logs/run_sim_parallel.log.

## What the full run records beyond the pilot
Per fit: REF and CAL (isotonic recalibration), the same fitted model evaluated on an
interpolation set inside the training window (interp_*), fill distance of the draw,
separation of the block from the training locations, covariate energy distance.
Per survey, species and edge: near oracle in sample Brier in the block and the extent
(the intrinsic term of the optimism decomposition) and a variogram practical range.
06b_theory_full.R turns these into the tables the manuscript's Section 4.4, 4.7 and
4.8 placeholders call for; 07_harmonise.R rebuilds the package at the new scale.

## Stages (run_full_local.R)
00 download (if needed) -> 01 preprocess (15 species/survey) -> 05 parallel simulation
-> 05b oracle -> 06 analysis -> 06b theory driven analyses -> 07 harmonisation and
package -> 08 figures -> 09 spatial figures.

## Reproducibility
A global seed drives a deterministic per cell seed, so every fit is reproducible
regardless of core count or order. Checkpoints are keyed to scale and replicate count
(results/parts_full_r100/), so pilot and full results never mix.
