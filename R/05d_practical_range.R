# =====================================================================
# 05d_practical_range.R
# Practical range of the residual variogram for every survey and species,
# profiled exactly over the range, with each species classed by whether its
# residual variogram identifies a range at all (R/variogram.R).
#
# The oracle stage (05b_oracle.R) fitted the same variogram with a single
# local optimiser start. For species whose residuals show no spatial
# structure the range was not identified and the optimiser returned its
# starting value; for species whose variogram has no sill within the survey
# the range ran away far beyond the survey extent. Either way the scaled
# separation of those species was an artefact. This stage repeats the fit
# with the same seed and the same subsample of hauls and reports, for each
# species, the profile least squares range and its regime.
#
#   results/practical_range.csv   dataset, species, range_km (identified
#     ranges only, NA otherwise), regime, sill_share, lag_min_km,
#     lag_max_km, dist_max_km, range_optim_km (the single start optimiser
#     result, which reproduces the practical_range_km column of
#     oracle_intrinsic.csv)
#
#   Rscript R/05d_practical_range.R
# =====================================================================
suppressMessages({ library(data.table); library(mgcv) })
if (!nzchar(Sys.getenv("SMD_ROOT"))) Sys.setenv(SMD_ROOT = getwd())
setwd(Sys.getenv("SMD_ROOT"))
source(file.path("R", "config.R")); source(file.path("R", "variogram.R"))

rows <- list()
for (code in names(DATASETS)) {
  r <- readRDS(file.path(PATHS$processed, paste0(code, "_processed.rds")))
  hauls <- as.data.table(r$hauls); env_used <- r$env_used
  for (si in seq_len(nrow(r$species))) {
    focal <- r$species$accepted_name[si]
    y <- r$occ[[focal]][match(hauls$haul_id, r$occ$haul_id)]
    full <- as.data.frame(hauls[, c("haul_id", "lon", "lat", "year", env_used), with = FALSE]); full$y <- y
    set.seed(GLOBAL_SEED + si)                     # as in 05b_oracle.R, so the same hauls are drawn
    v <- variogram_fit(full, env_used)
    rows[[length(rows) + 1]] <- data.table(dataset = code, species = focal, range_km = v$range_km, regime = v$regime,
                                           sill_share = v$sill_share, lag_min_km = v$lag_min_km,
                                           lag_max_km = v$lag_max_km, dist_max_km = v$dist_max_km, range_optim_km = v$range_optim_km)
  }
  cat(sprintf("[range] %s done\n", code))
}
pr <- rbindlist(rows)
fwrite(pr, file.path(PATHS$results, "practical_range.csv"))

# check against the unconstrained ranges recorded by the oracle stage
of <- file.path(PATHS$results, "oracle_intrinsic.csv")
if (file.exists(of)) {
  o <- unique(fread(of)[, .(dataset, species, oracle_km = practical_range_km)])
  ck <- merge(pr, o, by = c("dataset", "species"))
  cat(sprintf("[range] single start optimiser reproduces the oracle stage: max relative difference %.2e over %d species\n",
              ck[, max(abs(range_optim_km - oracle_km) / oracle_km, na.rm = TRUE)], nrow(ck)))
}
cat("[range] regimes:", paste(sprintf("%s %d", names(table(pr$regime)), table(pr$regime)), collapse = " | "), "\n")
cat("[range] written results/practical_range.csv\n")
