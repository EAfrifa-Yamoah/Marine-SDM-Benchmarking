# =====================================================================
# 06b_theory_full.R
# Analyses that need the quantities the full scale driver records; runs
# only when those columns are present (they are absent from the pilot).
#   A same fit optimism: one fitted predictor under both targets
#   B decomposition of optimism into intrinsic (oracle) and estimation terms
#   C Corollary 4: estimation optimism on scaled separation and covariate shift
#     (separation scaled by identified residual correlation ranges, 05d_practical_range.R)
#   D Proposition 5: fill distance in place of the coverage factor
#   E Corollary 5: Brier optimism on log n, same fit
#   F Theorem 2(ii): translation to the refinement term
#   G block separation against the residual correlation range (Corollaries 5 and 6)
# =====================================================================
suppressMessages({ library(data.table); library(lme4); library(mgcv) })
if (!nzchar(Sys.getenv("SMD_ROOT"))) Sys.setenv(SMD_ROOT = getwd()); setwd(Sys.getenv("SMD_ROOT"))
source(file.path("R", "config.R")); RES <- PATHS$results
say <- function(...) cat(sprintf(...), "\n")
fm_rds <- file.path(RES, "fit_metrics.rds")   # binary copy preferred (see 10_theory_analysis.R)
d <- (if (file.exists(fm_rds)) as.data.table(readRDS(fm_rds)) else fread(file.path(RES, "fit_metrics.csv")))[converged == TRUE]
if (!"interp_AUC" %in% names(d)) stop("fit_metrics.csv lacks the full scale columns; rerun the driver")
m <- d[target == "marginal" & is.finite(interp_AUC) & is.finite(AUC)]
m[, `:=`(opt_auc_samefit = interp_AUC - AUC, opt_brier_samefit = Brier - interp_Brier,
         cls = ifelse(spatial, "spatial", "non spatial"), log10n = log10(n))]

# A same fit optimism
A <- m[, .(opt_AUC = mean(opt_auc_samefit), opt_Brier = mean(opt_brier_samefit), n = .N), by = .(method, spatial)][order(-opt_AUC)]
A_ds <- m[, .(opt_AUC = mean(opt_auc_samefit), opt_Brier = mean(opt_brier_samefit)), by = .(dataset, cls)]
fwrite(A, file.path(RES, "table_samefit_optimism_by_method.csv")); fwrite(A_ds, file.path(RES, "table_samefit_optimism_by_dataset.csv"))
say("A same fit optimism (AUC): spatial %.3f | non spatial %.3f", m[spatial == TRUE, mean(opt_auc_samefit)], m[spatial == FALSE, mean(opt_auc_samefit)])

# B intrinsic and estimation terms (Brier scale)
of <- file.path(RES, "oracle_intrinsic.csv")
if (file.exists(of)) {
  o <- fread(of)
  # oracle proxy: median in sample Brier across methods in each region (flexible
  # learners over fit in sample, so the median is a conservative middle)
  ob <- o[, .(bs = median(brier_in_sample), range_optim_km = practical_range_km[1]), by = .(dataset, species, edge, region)]
  ow <- dcast(ob, dataset + species + edge + range_optim_km ~ region, value.var = "bs")
  ow[, intrinsic := block - extent]
  m <- merge(m, ow[, .(dataset, species, edge, intrinsic, range_optim_km)], by = c("dataset", "species", "edge"), all.x = TRUE)
  m[, estimation := opt_brier_samefit - intrinsic]
  B <- m[is.finite(intrinsic), .(intrinsic = mean(intrinsic), estimation_spatial = mean(estimation[spatial]),
                                  estimation_nonspatial = mean(estimation[!spatial]), opt_brier = mean(opt_brier_samefit)), by = dataset]
  fwrite(B, file.path(RES, "table_optimism_decomposition.csv")); print(B)
  say("B share of between survey spread in Brier optimism carried by the intrinsic term: %.2f", var(B$intrinsic) / var(B$opt_brier))
  # C Corollary 4: estimation term on separation in units of the residual correlation
  # range. The range comes from 05d_practical_range.R and is used only for the
  # species whose residual variogram identifies it; the original specification
  # (single start optimiser ranges for every species) and separation in
  # kilometres for every species are reported as sensitivity analyses.
  pr <- fread(file.path(RES, "practical_range.csv"))
  m <- merge(m, pr[, .(dataset, species, range_km, regime)], by = c("dataset", "species"), all.x = TRUE)
  m[, `:=`(sp_id = paste(dataset, species), cell_id = paste(dataset, species, n, coverage))]
  c4_fit <- function(dat, sep, spec) {
    dat <- copy(dat)[, sep_scaled := sep][is.finite(sep_scaled) & is.finite(cov_shift) & is.finite(estimation)]
    fit <- tryCatch(lmer(estimation ~ spatial * log1p(sep_scaled) + spatial * cov_shift + log10n + (1 | dataset) + (1 | sp_id) + (1 | cell_id),
                         data = dat, REML = TRUE), error = function(e) NULL)
    if (is.null(fit)) { say("C %s: model did not fit", spec); return(NULL) }
    cc <- summary(fit)$coefficients
    data.table(specification = spec, term = rownames(cc), estimate = cc[, 1], se = cc[, 2],
               n_fits = nrow(dat), n_species = uniqueN(dat$sp_id))
  }
  mi <- m[regime == "identified"]
  C <- c4_fit(mi, mi$separation_km / mi$range_km, "identified range")
  if (!is.null(C)) { fwrite(C[, .(term, estimate, se, n_fits, n_species)], file.path(RES, "table_corollary4_regression.csv")); print(C) }
  Cs <- rbind(c4_fit(m, m$separation_km / m$range_optim_km, "original: single start range, all species"),
              c4_fit(m, m$separation_km / 100, "separation per 100 km, all species"))
  fwrite(Cs, file.path(RES, "table_corollary4_sensitivity.csv")); print(Cs)
} else say("B/C skipped: oracle_intrinsic.csv not found (run 05b_oracle.R)")

# G separation of the held out block from its training locations, per marginal
# draw, against the residual correlation range (Corollaries 5 and 6, Section 4.8)
pr <- fread(file.path(RES, "practical_range.csv"))
dr <- unique(d[target == "marginal" & is.finite(separation_km), .(dataset, species, n, coverage, rep, separation_km)])
dr <- merge(dr, pr[, .(dataset, species, range_km, regime)], by = c("dataset", "species"))[, sp_id := paste(dataset, species)]
sep_summary <- function(x) x[, .(draws = .N, median_separation_km = median(separation_km),
  species_identified = uniqueN(sp_id[regime == "identified"]), draws_identified = sum(regime == "identified"),
  median_scaled_identified = median((separation_km / range_km)[regime == "identified"]),
  share_beyond_one_range_identified = mean((separation_km / range_km)[regime == "identified"] >= 1),
  species_no_structure = uniqueN(sp_id[regime == "no structure"]), species_below_lags = uniqueN(sp_id[regime == "below lags"]),
  species_no_sill = uniqueN(sp_id[regime == "no sill"]))]
G <- rbind(cbind(dataset = "all", sep_summary(dr)), dr[, sep_summary(.SD), by = dataset])
fwrite(G, file.path(RES, "table_block_separation.csv")); print(G)

# D fill distance in place of coverage (Proposition 5)
dm <- d[target == "marginal" & is.finite(fill_km) & fill_km > 0]
dm[, `:=`(log2_fill = log2(fill_km), sp_id = paste(dataset, species), cell_id = paste(dataset, species, n, coverage), log10n = log10(n))]
fit_d <- tryCatch(lmer(AUC ~ spatial + family + log10n + log2_fill + target_dummy + band + (1 | dataset) + (log2_fill | dataset) + (1 | sp_id) + (1 | cell_id),
              data = dm[, target_dummy := 0], REML = TRUE, control = lmerControl(check.conv.singular = .makeCC("ignore", tol = 1e-4))), error = function(e) NULL)
if (!is.null(fit_d)) {
fd <- summary(fit_d)$coefficients
say("D marginal AUC per doubling of fill distance: %.4f (se %.4f)", fd["log2_fill", 1], fd["log2_fill", 2])
fwrite(data.table(term = rownames(fd), estimate = fd[, 1], se = fd[, 2]), file.path(RES, "table_fill_distance_model.csv"))
fwrite(dm[, .(fill_km_mean = mean(fill_km), fill_km_sd = sd(fill_km)), by = .(dataset, coverage, n)], file.path(RES, "table_fill_distance_by_design.csv"))
} else say("D skipped: fill distance model did not fit on available data")

# E Corollary 5 on the same fit
fe <- lm(opt_brier_samefit ~ log10n * spatial, data = m); ce <- summary(fe)$coefficients
say("E slope of same fit Brier optimism on log10 n: non spatial %.4f | spatial %.4f | interaction p = %.3g",
    ce["log10n", 1], ce["log10n", 1] + ce["log10n:spatialTRUE", 1], ce["log10n:spatialTRUE", 4])
fwrite(data.table(term = rownames(ce), estimate = ce[, 1], se = ce[, 2], p = ce[, 4]), file.path(RES, "table_corollary5_samefit.csv"))

# F translation to the refinement term (Theorem 2(ii))
dr <- d[is.finite(AUC) & is.finite(REF) & is.finite(Brier)]
rng <- quantile(dr$AUC, c(.02, .98), names = FALSE); xg <- seq(rng[1], rng[2], length.out = 50)
g_ref <- gam(REF ~ s(AUC, k = 5), data = dr, method = "REML"); g_bs <- gam(Brier ~ s(AUC, k = 5), data = dr, method = "REML")
F <- data.table(target = c("REF", "Brier"),
                rel_resid_range = c(sd(residuals(g_ref)) / diff(range(dr$REF)), sd(residuals(g_bs)) / diff(range(dr$Brier))),
                rel_resid_idr = c(sd(residuals(g_ref)) / diff(quantile(dr$REF, c(.1, .9))), sd(residuals(g_bs)) / diff(quantile(dr$Brier, c(.1, .9)))))
fwrite(F, file.path(RES, "table_translation_to_refinement.csv")); print(F)
say("06b complete")
