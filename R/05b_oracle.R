# =====================================================================
# 05b_oracle.R
# Near oracle estimates for the decomposition of optimism (Theorem 4).
#
# For every survey, species and held out edge, each method is fitted to
# ALL available hauls in the held out block and, separately, to all hauls
# in the training extent; its in sample Brier score on each region
# approximates the irreducible uncertainty BS(p*, nu) there, up to the
# residual estimation error of a large sample fit. The difference between
# the two regions is the intrinsic term, common to every method; the
# estimation term of any fit follows by subtraction in 06_analysis.R.
#
# Also fitted once per survey and species: the practical range of an
# exponential variogram of the residuals of an environmental additive model,
# used to express block separation in correlation range units (Theorem 6).
#
# Fits are on large n and are capped at ORACLE_MAX_N hauls per region so
# that the stage stays affordable; that cap is recorded.
#   Rscript R/05b_oracle.R          (SMD_DATASETS subset honoured)
# =====================================================================
suppressMessages({ library(data.table); library(mgcv) })
if (!nzchar(Sys.getenv("SMD_ROOT"))) Sys.setenv(SMD_ROOT = getwd())
setwd(Sys.getenv("SMD_ROOT"))
for (f in c("config.R", "02_subsample.R", "03_methods.R", "04_metrics.R", "05_helpers.R"))
  source(file.path("R", f))
ORACLE_MAX_N <- as.integer(Sys.getenv("SMD_ORACLE_MAX_N", unset = "3000"))
ORACLE_METHODS <- setdiff(METHODS, "JointSDM")   # joint model needs the co species stack
sel <- Sys.getenv("SMD_DATASETS", unset = "")
ds_codes <- if (nzchar(sel)) trimws(strsplit(sel, ",")[[1]]) else names(DATASETS)
out_file <- file.path(PATHS$results, "oracle_intrinsic.csv")
done <- if (file.exists(out_file)) fread(out_file) else data.table()

# practical range of the residual variogram (km), bounded by the largest lag examined: R/variogram.R
source(file.path("R", "variogram.R"))

rows <- list()
for (code in ds_codes) {
  r <- readRDS(file.path(PATHS$processed, paste0(code, "_processed.rds")))
  hauls <- as.data.table(r$hauls); env_used <- r$env_used
  for (si in seq_len(nrow(r$species))) {
    focal <- r$species$accepted_name[si]
    y <- r$occ[[focal]][match(hauls$haul_id, r$occ$haul_id)]
    if (nrow(done) && any(done$dataset == code & done$species == focal)) next
    full <- as.data.frame(hauls[, c("haul_id", "lon", "lat", "year", env_used), with = FALSE]); full$y <- y
    set.seed(GLOBAL_SEED + si)
    rng <- practical_range(full, env_used)
    for (edge_id in 1:4) {
      hb <- .holdout_block(hauls$lon, hauls$lat, edge_id - 1L)   # rep_id %% 4 selects the edge
      for (region in c("block", "extent")) {
        idx <- if (region == "block") hb$test else hb$train_pool
        if (length(idx) > ORACLE_MAX_N) idx <- idx[sample.int(length(idx), ORACLE_MAX_N)]
        df <- full[idx, ]
        if (sum(df$y) < 5 || sum(1 - df$y) < 5) next
        for (mth in ORACLE_METHODS) {
          t0 <- Sys.time()
          res <- tryCatch(METHOD_FUNS[[mth]](df, df, env_used), error = function(e) NULL)
          if (is.null(res) || is.null(res$prob)) next
          rows[[length(rows) + 1]] <- data.table(
            dataset = code, species = focal, band = r$species$band[si], edge = hb$edge,
            region = region, method = mth, n_used = nrow(df), prevalence = mean(df$y),
            brier_in_sample = mean((df$y - res$prob)^2),
            practical_range_km = rng,
            fit_seconds = as.numeric(Sys.time() - t0, units = "secs"))
        }
      }
    }
    cat(sprintf("[oracle] %s %s done (%d rows so far)\n", code, focal, length(rows)))
    fwrite(rbindlist(c(list(done), rows), fill = TRUE), out_file)   # checkpoint per species
  }
}
cat("[oracle] complete:", out_file, "\n")
