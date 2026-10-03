# =====================================================================
# 06_analysis.R
# Analysis of the benchmark fit metrics.
#
# The pilot produces one row per fitted model. This script answers the
# five research questions of the study design from that table:
#
#   RQ1  how method performance scales with sample size across the data
#        limited range, and where methods break down
#   RQ2  the size of the spatial structure effect relative to the choice
#        of algorithm family
#   RQ3  the trade off between spatial coverage and total sample size
#   RQ4  validation optimism, the gap between conditional (interpolation)
#        and marginal (spatial extrapolation) performance
#   plus failure mode characterisation: split failures, miscalibration
#        and overfitting
#
# Hierarchical linear models separate design effects from the nested
# structure of replicates within design cells within species within
# datasets. Outputs are written as tidy CSV tables and a single analysis
# object consumed by the figure and reporting scripts.
# =====================================================================

suppressMessages({ library(data.table); library(lme4) })
source(file.path("R", "config.R"))

RES <- PATHS$results
d <- as.data.table(readRDS(file.path(RES, "fit_metrics.rds")))
fails <- if (file.exists(file.path(RES, "split_failures.csv")))
  fread(file.path(RES, "split_failures.csv")) else data.table()

# derived analysis variables
d[, log10n := log10(n)]
d[, cell_id := paste(dataset, species, n, coverage, target, sep = "|")]
d[, sp_id := paste(dataset, species, sep = "|")]
d[, method := factor(method, levels = METHODS)]
d[, band := factor(band, levels = c("rare", "intermediate", "common"))]
d[, coverage := factor(coverage,
      levels = c("clustered", "intermediate", "distributed"))]
# absolute miscalibration on the log slope scale (0 = perfect slope of 1)
d[, miscal := abs(log(pmax(calib_slope, 1e-3)))]

analysis <- list()
analysis$n_fits <- nrow(d)
analysis$methods <- METHODS
analysis$converged_rate <- mean(d$converged)

log_file <- file.path(PATHS$logs, "analysis.log")
say <- function(...) {
  msg <- sprintf(...); cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}
cat("", file = log_file)
say("=== analysis on %d fits, %d methods ===", nrow(d), length(METHODS))

# ---------------------------------------------------------------------
# Overall performance table, by method and inference target
# ---------------------------------------------------------------------
perf <- d[, .(
  AUC = mean(AUC, na.rm = TRUE), AUC_sd = sd(AUC, na.rm = TRUE),
  TSS = mean(TSS, na.rm = TRUE), RMSE = mean(RMSE, na.rm = TRUE),
  R2 = mean(R2, na.rm = TRUE), Brier = mean(Brier, na.rm = TRUE),
  Spearman = mean(Spearman, na.rm = TRUE),
  calib_slope = median(calib_slope, na.rm = TRUE),
  gap_AUC = mean(gap_AUC, na.rm = TRUE), n_fits = .N),
  by = .(method, target)]
setorder(perf, target, -AUC)
fwrite(perf, file.path(RES, "table_performance_by_method_target.csv"))
analysis$perf <- perf

# ---------------------------------------------------------------------
# RQ1  sample size scaling
# per method AUC against log10 sample size, slope and marginal means
# ---------------------------------------------------------------------
scaling <- d[, .(AUC = mean(AUC, na.rm = TRUE), sd = sd(AUC, na.rm = TRUE),
                 n_fits = .N), by = .(method, n)]
setorder(scaling, method, n)
fwrite(scaling, file.path(RES, "table_scaling_by_method_n.csv"))
analysis$scaling <- scaling

# per method slope of AUC on log10(n), conditional target only, so the
# scaling is read on the interpolation performance the design defines it on
slopes <- d[target == "conditional",
  {
    fit <- lm(AUC ~ log10n)
    .(slope_per_decade = unname(coef(fit)[2]),
      auc_at_n30  = predict(fit, data.frame(log10n = log10(30))),
      auc_at_n500 = predict(fit, data.frame(log10n = log10(500))))
  }, by = method]
setorder(slopes, -slope_per_decade)
fwrite(slopes, file.path(RES, "table_scaling_slopes.csv"))
analysis$scaling_slopes <- slopes
say("RQ1 sample size: pooled AUC %.3f (n=30) to %.3f (n=500)",
    scaling[n == 30, mean(AUC)], scaling[n == 500, mean(AUC)])

# ---------------------------------------------------------------------
# RQ2  spatial structure effect versus algorithm family
# hierarchical model with a spatial indicator and method family, random
# intercepts for the nested design. Variance partition quantifies how
# much of the spread each source explains.
# ---------------------------------------------------------------------
d[, family := factor(family)]
# A random intercept for dataset shifts the mean performance of a survey but
# cannot represent a survey specific validation optimism. The dataset by
# target term supplies exactly that, and makes the heterogeneity reported in
# Section 3.7 an estimable variance component rather than an invisible one.
m_main <- lmer(
  AUC ~ spatial + family + log10n + coverage + target + band +
    (1 | dataset) + (1 | dataset:target) + (1 | sp_id) + (1 | cell_id),
  data = d, REML = TRUE,
  control = lmerControl(optimizer = "bobyqa",
                        check.conv.singular = .makeCC("ignore", tol = 1e-4)))
vc <- as.data.frame(VarCorr(m_main))
vc <- vc[, c("grp", "vcov")]
vc$prop <- vc$vcov / sum(vc$vcov)
fwrite(as.data.table(vc), file.path(RES, "table_variance_partition.csv"))
analysis$variance_partition <- vc

fe <- summary(m_main)$coefficients
analysis$fixed_effects <- as.data.table(fe, keep.rownames = "term")
fwrite(analysis$fixed_effects, file.path(RES, "table_fixed_effects.csv"))

# spatial effect: contrast of spatial versus non spatial methods within
# each target, on AUC, as a directly interpretable mean difference
sp_eff <- d[, .(AUC = mean(AUC, na.rm = TRUE)), by = .(spatial, target)]
sp_eff <- dcast(sp_eff, target ~ spatial, value.var = "AUC")
setnames(sp_eff, c("FALSE", "TRUE"), c("non_spatial", "spatial"))
sp_eff[, spatial_advantage := spatial - non_spatial]
fwrite(sp_eff, file.path(RES, "table_spatial_effect.csv"))
analysis$spatial_effect <- sp_eff
say("RQ2 spatial advantage: conditional %+.3f AUC, marginal %+.3f AUC",
    sp_eff[target == "conditional", spatial_advantage],
    sp_eff[target == "marginal", spatial_advantage])
say("RQ2 variance explained: cell %.0f%%, species %.0f%%, dataset %.0f%%, dataset:target %.0f%%, residual %.0f%%",
    100 * vc$prop[vc$grp == "cell_id"], 100 * vc$prop[vc$grp == "sp_id"],
    100 * vc$prop[vc$grp == "dataset"],
    100 * sum(vc$prop[vc$grp == "dataset:target"]),
    100 * vc$prop[vc$grp == "Residual"])

# ---------------------------------------------------------------------
# RQ3  coverage versus sample size trade off
# mean AUC over the coverage by sample size grid, and the sample size
# equivalent of moving from clustered to distributed coverage
# ---------------------------------------------------------------------
cov_n <- d[target == "marginal",
  .(AUC = mean(AUC, na.rm = TRUE), n_fits = .N), by = .(coverage, n)]
setorder(cov_n, coverage, n)
cov_grid <- dcast(cov_n, n ~ coverage, value.var = "AUC")
fwrite(cov_grid, file.path(RES, "table_coverage_by_n.csv"))
analysis$coverage_grid <- cov_grid

# coverage main effect at marginal target (extrapolation is where spatial
# coverage should matter most)
cov_eff <- d[target == "marginal",
             .(AUC = mean(AUC, na.rm = TRUE)), by = coverage]
analysis$coverage_effect <- cov_eff
say("RQ3 coverage at marginal: clustered %.3f, intermediate %.3f, distributed %.3f",
    cov_eff[coverage == "clustered", AUC],
    cov_eff[coverage == "intermediate", AUC],
    cov_eff[coverage == "distributed", AUC])

# ---------------------------------------------------------------------
# RQ4  validation optimism, conditional minus marginal
# paired within the same dataset, species, sample size, coverage, replicate
# and method so the contrast isolates the inference target
# ---------------------------------------------------------------------
key <- c("dataset", "sp_id", "species", "band", "n", "coverage",
         "rep", "method", "family", "spatial")
wide <- dcast(d, dataset + sp_id + species + band + n + coverage + rep +
                method + family + spatial ~ target,
              value.var = "AUC")
wide <- wide[!is.na(conditional) & !is.na(marginal)]
wide[, optimism := conditional - marginal]
opt_method <- wide[, .(optimism = mean(optimism),
                       opt_sd = sd(optimism),
                       cond = mean(conditional), marg = mean(marginal),
                       n_pairs = .N), by = .(method, spatial)]
setorder(opt_method, -optimism)
fwrite(opt_method, file.path(RES, "table_validation_optimism.csv"))
analysis$optimism <- opt_method

# optimism by sample size: does more data reduce the optimism gap
opt_n <- wide[, .(optimism = mean(optimism)), by = n][order(n)]
analysis$optimism_by_n <- opt_n

# optimism by coverage: does broader coverage reduce extrapolation penalty
opt_cov <- wide[, .(optimism = mean(optimism)), by = coverage]
analysis$optimism_by_coverage <- opt_cov
say("RQ4 optimism: spatial methods %+.3f AUC, non spatial %+.3f AUC",
    wide[spatial == TRUE, mean(optimism)],
    wide[spatial == FALSE, mean(optimism)])

# ---------------------------------------------------------------------
# Heterogeneity across surveys
# The pooled figures above average over four surveys whose sampling
# profiles and geography differ (Section 2.3). Here every headline
# contrast is recomputed within survey, so the reader can see how much
# the sampling profile conditions the result rather than inferring it
# from a single dataset variance component.
# ---------------------------------------------------------------------

# (a) validation optimism by survey, paired, split by method class
opt_ds <- wide[, .(conditional = mean(conditional),
                   marginal    = mean(marginal),
                   optimism    = mean(optimism),
                   opt_se      = sd(optimism) / sqrt(.N),
                   n_pairs     = .N), by = .(dataset, spatial)]
opt_ds_all <- wide[, .(conditional = mean(conditional),
                       marginal    = mean(marginal),
                       optimism    = mean(optimism),
                       opt_se      = sd(optimism) / sqrt(.N),
                       n_pairs     = .N), by = .(dataset)]
opt_ds_all[, spatial := NA]
setorder(opt_ds_all, -optimism)
fwrite(opt_ds,     file.path(RES, "table_optimism_by_dataset.csv"))
fwrite(opt_ds_all, file.path(RES, "table_optimism_by_dataset_pooled.csv"))
analysis$optimism_by_dataset <- opt_ds
analysis$optimism_by_dataset_pooled <- opt_ds_all

# (b) spatial advantage by survey and target: does the sign change
sp_ds <- d[, .(AUC = mean(AUC, na.rm = TRUE)), by = .(dataset, target, spatial)]
sp_ds <- dcast(sp_ds, dataset + target ~ spatial, value.var = "AUC")
setnames(sp_ds, c("FALSE", "TRUE"), c("non_spatial", "spatial"))
sp_ds[, spatial_advantage := spatial - non_spatial]
setorder(sp_ds, target, -spatial_advantage)
fwrite(sp_ds, file.path(RES, "table_spatial_effect_by_dataset.csv"))
analysis$spatial_effect_by_dataset <- sp_ds

# (c) coverage lift by survey, under the marginal target
cov_ds <- d[target == "marginal",
            .(AUC = mean(AUC, na.rm = TRUE)), by = .(dataset, coverage)]
cov_ds <- dcast(cov_ds, dataset ~ coverage, value.var = "AUC")
setcolorder(cov_ds, c("dataset", "clustered", "intermediate", "distributed"))
cov_ds[, coverage_lift := distributed - clustered]
setorder(cov_ds, -coverage_lift)
fwrite(cov_ds, file.path(RES, "table_coverage_by_dataset.csv"))
analysis$coverage_by_dataset <- cov_ds

# (d) survey sampling profile descriptors, tying the numbers above back to
# the point clouds in Figure 1. Station revisitation is reported because
# the surveys pool many years, so a nearest neighbour statistic computed on
# the pooled cloud measures repeat visits rather than the spatial design.
prof <- rbindlist(lapply(names(DATASETS), function(code) {
  pf <- file.path(PATHS$processed, paste0(code, "_processed.rds"))
  if (!file.exists(pf)) return(NULL)
  h <- as.data.table(readRDS(pf)$hauls)
  latm <- mean(h$lat)
  h[, `:=`(sx = round(lon, 2), sy = round(lat, 2))]
  u <- unique(h, by = c("sx", "sy"))
  cs <- 0.25
  cells <- unique(data.table(cx = floor(h$lon / cs), cy = floor(h$lat / cs)))
  cell_km2 <- (cs * 111 * cos(latm * pi / 180)) * (cs * 111)
  data.table(
    dataset          = code,
    hauls            = nrow(h),
    years            = length(unique(h$year)),
    unique_stations  = nrow(u),
    hauls_per_station= round(nrow(h) / nrow(u), 2),
    extent_km2       = round(nrow(cells) * cell_km2),
    span_x_km        = round(diff(range(h$lon)) * 111 * cos(latm * pi / 180)),
    span_y_km        = round(diff(range(h$lat)) * 111),
    n_env            = length(readRDS(pf)$env_used),
    env              = paste(readRDS(pf)$env_used, collapse = "+"))
}))
fwrite(prof, file.path(RES, "table_survey_profiles.csv"))
analysis$survey_profiles <- prof

say("HETEROGENEITY optimism by survey: %s",
    paste(sprintf("%s %.3f", opt_ds_all$dataset, opt_ds_all$optimism),
          collapse = " | "))
say("HETEROGENEITY spatial advantage sign flips: %d of %d survey by target cells negative",
    sp_ds[spatial_advantage < 0, .N], nrow(sp_ds))

# ---------------------------------------------------------------------
# Failure modes
# split failures (degenerate occurrence), miscalibration, overfitting
# ---------------------------------------------------------------------
if (nrow(fails)) {
  fail_tab <- fails[, .(n_fail = .N),
                    by = .(band, n, coverage)][order(band, n)]
  fail_by_band <- fails[, .(n_fail = .N), by = band]
  # attempted split cells per band for a rate: cells x replicates
  attempted <- d[, .(fits = .N), by = .(band, n, coverage)]
  fwrite(fail_tab, file.path(RES, "table_split_failures.csv"))
  analysis$split_failures <- fail_tab
  analysis$split_failures_by_band <- fail_by_band
  say("failures: %d degenerate splits, %s",
      nrow(fails),
      paste(sprintf("%s=%d", fail_by_band$band, fail_by_band$n_fail),
            collapse = ", "))
} else analysis$split_failures <- data.table()

# miscalibration by method: median calibration slope and share of fits
# whose slope departs strongly from one (over or under confidence)
miscal <- d[, .(
  calib_slope_med = median(calib_slope, na.rm = TRUE),
  frac_overconfident = mean(calib_slope < 0.8, na.rm = TRUE),
  frac_underconfident = mean(calib_slope > 1.25, na.rm = TRUE),
  miscal_med = median(miscal, na.rm = TRUE)),
  by = method]
setorder(miscal, -miscal_med)
fwrite(miscal, file.path(RES, "table_calibration_by_method.csv"))
analysis$calibration <- miscal

# overfitting: train minus test AUC gap by method and sample size
overfit <- d[, .(gap_AUC = mean(gap_AUC, na.rm = TRUE)),
             by = .(method, n)]
overfit_method <- d[, .(gap_AUC = mean(gap_AUC, na.rm = TRUE),
                        gap_sd = sd(gap_AUC, na.rm = TRUE)), by = method]
setorder(overfit_method, -gap_AUC)
fwrite(overfit_method, file.path(RES, "table_overfitting_by_method.csv"))
analysis$overfitting <- overfit
analysis$overfitting_method <- overfit_method
say("overfitting: largest train-test AUC gap %s (%.3f), smallest %s (%.3f)",
    overfit_method[1, method], overfit_method[1, gap_AUC],
    overfit_method[.N, method], overfit_method[.N, gap_AUC])

saveRDS(analysis, file.path(RES, "analysis.rds"))
say("analysis object written to results/analysis.rds")
