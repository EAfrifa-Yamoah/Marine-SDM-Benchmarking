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

Observed on a 2 core sandbox (R 4.3, mgcv 1.9.3): 0.12 s per fit on average, 1,177,771
scored fits in about 36 core hours; oracle stage 65 min; 06_analysis 12 min (lme4 on
1.18 million rows, 2 GB); 07_harmonise 70 min (210 REML smooths); 08 and 09 figures
10 min each; 10_theory_analysis 50 min before the bootstrap, with a peak of 6 GB when
another R process runs alongside it, so run it alone on an 8 GB machine. Its cluster
bootstrap is about 13 min per replicate at this scale (100 replicates, resumable within
a replicate; set SMD_BOOT_SECS to the wall clock you can give one invocation and re run
until "bootstrap replicates completed: 100 of 100").

When the parallel driver cannot be used (no sockets, or a time boxed job queue), run
R/05_run_sim_shard.R with SMD_SHARDS and SMD_SHARD for each shard and SMD_TIME_BUDGET in
seconds; every replicate block is checkpointed, and R/05c_merge_parts.R merges the
checkpoints into results/fit_metrics.rds.

## What the full run records beyond the pilot
Per fit: REF and CAL (isotonic recalibration), the same fitted model evaluated on an
interpolation set inside the training window (interp_*), fill distance of the draw,
separation of the block from the training locations, covariate energy distance.
Per survey, species and edge: near oracle in sample Brier in the block and the extent
(the intrinsic term of the optimism decomposition). 05d_practical_range.R fits the residual variogram of
every species by least squares profiled exactly over the range and records whether the range is
identified; only identified ranges scale block separation (results/practical_range.csv).
06b_theory_full.R turns these into the tables the manuscript's Section 4.4, 4.7 and
4.8 placeholders call for; 07_harmonise.R rebuilds the package at the new scale.

## Stages (run_full_local.R)
00 download (if needed) -> 01 preprocess (15 species/survey) -> 05 parallel simulation
-> 05b oracle -> 05d practical range -> 06 analysis -> 06b theory driven analyses -> 07 harmonisation and
package -> 08 figures -> 09 spatial figures.

After run_full_local.R, with SMD_SCALE=full and SMD_REPLICATES=100:
10 theory analysis and species bootstrap (repeat until 100 of 100 replicates) -> 07b refresh
the package's theory table -> 11 manuscript numbers (results/manuscript_numbers.csv) ->
12 Supporting Information (results/supplement/, figures/supplement/). 12 takes under a
minute; both read fit_metrics.rds, not the csv, which network backed storage can leave stale.

## Reproducibility
A global seed drives a deterministic per cell seed, so every fit is reproducible
regardless of core count or order. Checkpoints are keyed to scale and replicate count
(results/parts_full_r100/), so pilot and full results never mix.
