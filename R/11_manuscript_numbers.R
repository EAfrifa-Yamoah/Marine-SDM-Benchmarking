# =====================================================================
# 11_manuscript_numbers.R
# Every number the manuscript cites, recomputed from the pipeline outputs,
# with its definition and the file it comes from. Run after 06 to 10.
#
#   results/manuscript_numbers.csv   key, value, shown, definition, source
#
# The column `shown` is the value exactly as it is printed in the
# manuscript. After reproducing the run locally, compare this file with
# the highlighted numbers in the manuscript; any mismatch is a number that
# must not go into the paper.
#
# Two quantities are computed here rather than read from a table:
#   (1) contrasts that exclude JointSDM, reported as a sensitivity analysis
#       because the joint model analogue shares one environmental response
#       across species (see R/03_methods.R);
#   (2) standard errors clustered by species for the sample size slopes of
#       Corollary 5, because the fits within a species are not independent
#       and ordinary least squares standard errors overstate precision.
# It also rewrites table_spearman_tie_effect_by_method.csv, whose per method
# counts were not split by method in earlier versions of 10_theory_analysis.R.
# =====================================================================
suppressMessages({ library(data.table); library(sandwich) })
if (!nzchar(Sys.getenv("SMD_ROOT"))) Sys.setenv(SMD_ROOT = getwd())
setwd(Sys.getenv("SMD_ROOT")); source(file.path("R", "config.R"))
RES <- PATHS$results
tab <- function(f) fread(file.path(RES, f))

N <- list()
add <- function(key, value, shown, definition, source) {
  N[[length(N) + 1]] <<- data.table(id = key, value = as.numeric(value),
                                    shown = shown, definition = definition,
                                    source = source)
  invisible(value)
}
f3 <- function(x) sprintf("%.3f", x)
fs <- function(x) sprintf("%+.3f", x)                 # signed, three decimals
f4 <- function(x) sprintf("%.4f", x)
fc <- function(x) format(round(x), big.mark = ",", scientific = FALSE)
fp <- function(x, d = 1) sprintf(paste0("%.", d, "f%%"), 100 * x)

d0 <- as.data.table(readRDS(file.path(RES, "fit_metrics.rds")))
d0[, sp_id := paste(dataset, species)]
dc <- d0[converged == TRUE & is.finite(AUC) & is.finite(R2) & is.finite(Brier)]
SM <- c("SpatialRF", "SpatialGAM", "GeostatGP")      # spatial methods other than JointSDM

# ---------------------------------------------------------------------
# design and execution
# ---------------------------------------------------------------------
sf <- tab("split_failures.csv")
planned <- length(DATASETS) * N_SPECIES_PER_DATASET * 30 * N_REPLICATES * length(METHODS)
add("fits_planned", planned, fc(planned), "surveys x species x design cells x replicates x methods", "config.R")
add("splits_degenerate", nrow(sf), fc(nrow(sf)), "replicate splits with a degenerate occurrence vector", "split_failures.csv")
add("fits_lost", nrow(sf) * length(METHODS), fc(nrow(sf) * length(METHODS)), "fits lost to degenerate splits", "split_failures.csv")
add("fits_scored", nrow(d0), fc(nrow(d0)), "scored fits", "fit_metrics.rds")
add("fits_converged", nrow(dc), fc(nrow(dc)), "converged fits with finite AUC, R2 and Brier", "fit_metrics.rds")
add("conv_rate", mean(d0$converged), fp(mean(d0$converged), 2), "share of scored fits that converged", "fit_metrics.rds")
for (b in c("rare", "intermediate", "common"))
  add(paste0("degenerate_", b), sf[band == b, .N], fc(sf[band == b, .N]), paste("degenerate splits,", b, "band"), "split_failures.csv")
for (nn in SAMPLE_SIZES)
  add(paste0("degenerate_n", nn), sf[n == nn, .N], fc(sf[n == nn, .N]), paste("degenerate splits at n =", nn), "split_failures.csv")
add("splits_attempted", planned / length(METHODS), fc(planned / length(METHODS)), "replicate splits attempted", "config.R")
add("degenerate_share_rare", sf[band == "rare", .N] / nrow(sf), fp(sf[band == "rare", .N] / nrow(sf), 0), "share of degenerate splits in the rare band", "split_failures.csv")
add("degenerate_share_n50",sf[n <= 50, .N] / nrow(sf), fp(sf[n <= 50, .N] / nrow(sf), 0), "share of degenerate splits at n = 30 or 50", "split_failures.csv")
add("species_total", uniqueN(d0$sp_id), fc(uniqueN(d0$sp_id)), "survey by species combinations", "fit_metrics.rds")
add("test_size", median(d0$n_test), fc(median(d0$n_test)), "test hauls per fit (median)", "fit_metrics.rds")
add("fit_seconds_mean", mean(d0$fit_seconds, na.rm = TRUE), sprintf("%.2f", mean(d0$fit_seconds, na.rm = TRUE)), "mean fit time, seconds", "fit_metrics.rds")
add("core_hours", sum(d0$fit_seconds, na.rm = TRUE) / 3600, fc(sum(d0$fit_seconds, na.rm = TRUE) / 3600), "summed fit time, hours", "fit_metrics.rds")
add("oracle_fits", nrow(tab("oracle_intrinsic.csv")), fc(nrow(tab("oracle_intrinsic.csv"))), "near oracle fits", "oracle_intrinsic.csv")
prev_rng <- d0[, .(p = prevalence[1]), by = .(sp_id, band)]
for (b in c("rare", "intermediate", "common")) {
  add(paste0("prev_min_", b), prev_rng[band == b, min(p)], fp(prev_rng[band == b, min(p)]), paste("lowest survey prevalence,", b, "band"), "fit_metrics.rds")
  add(paste0("prev_max_", b), prev_rng[band == b, max(p)], fp(prev_rng[band == b, max(p)]), paste("highest survey prevalence,", b, "band"), "fit_metrics.rds")
}

# ---------------------------------------------------------------------
# RQ1 sample size
# ---------------------------------------------------------------------
for (tg in c("conditional", "marginal")) for (nn in c(30, 500)) {
  v <- d0[target == tg & n == nn, mean(AUC, na.rm = TRUE)]
  add(sprintf("auc_%s_n%d", tg, nn), v, f3(v), sprintf("mean test AUC, %s target, n = %d, all methods", tg, nn), "fit_metrics.rds")
}
cm <- d0[target == "conditional", .(AUC = mean(AUC, na.rm = TRUE)), by = .(method, n)]
for (m in METHODS) for (nn in c(30, 500))
  add(sprintf("auc_cond_%s_n%d", m, nn), cm[method == m & n == nn, AUC], f3(cm[method == m & n == nn, AUC]),
      sprintf("mean conditional AUC, %s, n = %d", m, nn), "fit_metrics.rds")
sl <- tab("table_scaling_slopes.csv")
for (m in METHODS) add(paste0("slope_", m), sl[method == m, slope_per_decade], f3(sl[method == m, slope_per_decade]),
                       paste("conditional AUC gain per decade of n,", m), "table_scaling_slopes.csv")
cross <- function(m) { w <- dcast(cm[method %in% c(m, "RF")], n ~ method, value.var = "AUC"); w[get(m) > RF, min(n)] }
for (m in c("SpatialGAM", "GeostatGP")) add(paste0("crossover_", m), cross(m), fc(cross(m)),
                                            paste("smallest n at which", m, "exceeds RF under the conditional target"), "fit_metrics.rds")
pf <- tab("table_performance_by_method_target.csv")
for (m in METHODS) for (tg in c("conditional", "marginal"))
  add(sprintf("auc_%s_%s", tg, m), pf[method == m & target == tg, AUC], f3(pf[method == m & target == tg, AUC]),
      sprintf("mean AUC, %s target, %s", tg, m), "table_performance_by_method_target.csv")

# ---------------------------------------------------------------------
# RQ2 spatial structure, variance partition
# ---------------------------------------------------------------------
se <- tab("table_spatial_effect.csv")
add("spadv_cond", se[target == "conditional", spatial_advantage], fs(se[target == "conditional", spatial_advantage]), "spatial minus non spatial mean AUC, conditional", "table_spatial_effect.csv")
add("spadv_marg", se[target == "marginal", spatial_advantage], fs(se[target == "marginal", spatial_advantage]), "spatial minus non spatial mean AUC, marginal", "table_spatial_effect.csv")
for (tg in c("conditional", "marginal")) {
  v <- d0[target == tg & method %in% SM, mean(AUC, na.rm = TRUE)] - d0[target == tg & spatial == FALSE, mean(AUC, na.rm = TRUE)]
  add(paste0("spadv_exJ_", substr(tg, 1, 4)), v, fs(v), paste("SpatialRF, SpatialGAM, GeostatGP minus non spatial mean AUC,", tg), "fit_metrics.rds")
}
fe <- tab("table_fixed_effects.csv")
add("fe_spatial", fe[term == "spatialTRUE", Estimate], fs(fe[term == "spatialTRUE", Estimate]), "spatial fixed effect, AUC model adjusted for family", "table_fixed_effects.csv")
add("fe_spatial_se", fe[term == "spatialTRUE", `Std. Error`], f4(fe[term == "spatialTRUE", `Std. Error`]), "its standard error", "table_fixed_effects.csv")
add("fe_jsdm", fe[term == "familyjsdm", Estimate], fs(fe[term == "familyjsdm", Estimate]), "joint model family effect", "table_fixed_effects.csv")
add("fe_cov_int", fe[term == "coverageintermediate", Estimate], fs(fe[term == "coverageintermediate", Estimate]), "intermediate coverage effect", "table_fixed_effects.csv")
add("fe_cov_dist", fe[term == "coveragedistributed", Estimate], fs(fe[term == "coveragedistributed", Estimate]), "distributed coverage effect", "table_fixed_effects.csv")
add("fe_log10n", fe[term == "log10n", Estimate], f3(fe[term == "log10n", Estimate]), "AUC per decade of n, adjusted", "table_fixed_effects.csv")
add("fe_marginal", fe[term == "targetmarginal", Estimate], fs(fe[term == "targetmarginal", Estimate]), "marginal target effect", "table_fixed_effects.csv")
vp <- tab("table_variance_partition.csv")
lab <- c(Residual = "residual", sp_id = "species", dataset = "survey", cell_id = "cell", `dataset:target` = "survey_target")
for (g in names(lab)) add(paste0("vp_", lab[g]), vp[grp == g, prop], fp(vp[grp == g, prop], 0), paste("variance share,", lab[g]), "table_variance_partition.csv")
vs <- vp[grp %in% c("dataset", "dataset:target"), sum(prop)]
add("vp_survey_total", vs, fp(vs, 0), "survey and survey by target shares combined", "table_variance_partition.csv")
add("vp_species_floor", vp[grp == "sp_id", prop] / N_SPECIES_PER_DATASET, fp(vp[grp == "sp_id", prop] / N_SPECIES_PER_DATASET, 1),
    "species share divided by species per survey (noise floor of a survey mean)", "table_variance_partition.csv")

# ---------------------------------------------------------------------
# RQ3 coverage
# ---------------------------------------------------------------------
cv <- d0[target == "marginal", .(AUC = mean(AUC, na.rm = TRUE)), by = coverage]
for (c in c("clustered", "intermediate", "distributed")) add(paste0("cov_", c), cv[coverage == c, AUC], f3(cv[coverage == c, AUC]), paste("mean marginal AUC,", c), "fit_metrics.rds")
add("cov_lift", cv[coverage == "distributed", AUC] - cv[coverage == "clustered", AUC], f3(cv[coverage == "distributed", AUC] - cv[coverage == "clustered", AUC]), "distributed minus clustered marginal AUC", "fit_metrics.rds")
cg <- tab("table_coverage_by_n.csv")
for (nn in SAMPLE_SIZES) for (c in c("clustered", "intermediate", "distributed"))
  add(sprintf("covgrid_%s_n%d", c, nn), cg[n == nn, get(c)], f3(cg[n == nn, get(c)]), sprintf("mean marginal AUC, %s, n = %d", c, nn), "table_coverage_by_n.csv")

# ---------------------------------------------------------------------
# RQ4 validation optimism
# ---------------------------------------------------------------------
w <- dcast(d0, dataset + sp_id + band + n + coverage + rep + method + spatial ~ target, value.var = "AUC")
w <- w[!is.na(conditional) & !is.na(marginal)][, opt := conditional - marginal]
vo <- tab("table_validation_optimism.csv")
for (m in METHODS) {
  add(paste0("opt_", m), vo[method == m, optimism], f3(vo[method == m, optimism]), paste("conditional minus marginal AUC, paired,", m), "table_validation_optimism.csv")
  add(paste0("opt_se_", m), vo[method == m, opt_sd / sqrt(n_pairs)], f4(vo[method == m, opt_sd / sqrt(n_pairs)]), paste("its standard error,", m), "table_validation_optimism.csv")
}
# standard errors from species means (fits within a species are not independent)
spm <- w[, .(o = mean(opt)), by = .(method, sp_id)][, .(se = sd(o) / sqrt(.N)), by = method]
add("opt_se_cl_min", min(spm$se), f3(min(spm$se)), "smallest species clustered standard error of a method's optimism", "fit_metrics.rds")
add("opt_se_cl_max", max(spm$se), f3(max(spm$se)), "largest species clustered standard error of a method's optimism", "fit_metrics.rds")
sps <- w[, .(o = mean(opt)), by = .(dataset, sp_id)][, .(se = sd(o) / sqrt(.N)), by = dataset]
for (s in sps$dataset) add(paste0("optds_se_cl_", gsub("-", "", s)), sps[dataset == s, se], f3(sps[dataset == s, se]),
                           paste("species clustered standard error of survey optimism,", s), "fit_metrics.rds")
add("optds_se_cl_min", min(sps$se), f3(min(sps$se)), "smallest species clustered standard error of a survey's optimism", "fit_metrics.rds")
add("optds_se_cl_max", max(sps$se), f3(max(sps$se)), "largest species clustered standard error of a survey's optimism", "fit_metrics.rds")
hr <- tab("table_harmonisation_reliability.csv")
add("spearman_auc_r2", hr[from == "AUC" & to == "R2", spearman], sprintf("%.2f", hr[from == "AUC" & to == "R2", spearman]), "Spearman correlation between AUC and R2 over fits", "table_harmonisation_reliability.csv")
add("opt_spatial", w[spatial == TRUE, mean(opt)], f3(w[spatial == TRUE, mean(opt)]), "optimism, spatial class including JointSDM", "fit_metrics.rds")
add("opt_nonspatial", w[spatial == FALSE, mean(opt)], f3(w[spatial == FALSE, mean(opt)]), "optimism, non spatial class", "fit_metrics.rds")
add("opt_spatial_exJ", w[method %in% SM, mean(opt)], f3(w[method %in% SM, mean(opt)]), "optimism, SpatialRF, SpatialGAM and GeostatGP", "fit_metrics.rds")
add("opt_excess_exJ", w[method %in% SM, mean(opt)] - w[spatial == FALSE, mean(opt)], f3(w[method %in% SM, mean(opt)] - w[spatial == FALSE, mean(opt)]), "excess optimism of the three spatial smoothers over non spatial methods", "fit_metrics.rds")
on <- w[, .(o = mean(opt)), by = n][order(n)]
for (nn in SAMPLE_SIZES) add(paste0("opt_n", nn), on[n == nn, o], f3(on[n == nn, o]), paste("pooled optimism at n =", nn), "fit_metrics.rds")
oc <- w[, .(o = mean(opt)), by = coverage]
for (c in c("clustered", "intermediate", "distributed")) add(paste0("opt_cov_", c), oc[coverage == c, o], f3(oc[coverage == c, o]), paste("pooled optimism,", c, "coverage"), "fit_metrics.rds")
mm <- d0[converged == TRUE & target == "marginal" & is.finite(interp_AUC) & is.finite(AUC)]
mm[, sf_opt := interp_AUC - AUC]
add("sfopt_spatial", mm[spatial == TRUE, mean(sf_opt)], f3(mm[spatial == TRUE, mean(sf_opt)]), "same fit optimism, spatial class", "fit_metrics.rds")
add("sfopt_nonspatial", mm[spatial == FALSE, mean(sf_opt)], f3(mm[spatial == FALSE, mean(sf_opt)]), "same fit optimism, non spatial class", "fit_metrics.rds")
add("sfopt_spatial_exJ", mm[method %in% SM, mean(sf_opt)], f3(mm[method %in% SM, mean(sf_opt)]), "same fit optimism, three spatial smoothers", "fit_metrics.rds")
sfm <- tab("table_samefit_optimism_by_method.csv")
for (m in METHODS) add(paste0("sfopt_", m), sfm[method == m, opt_AUC], f3(sfm[method == m, opt_AUC]), paste("same fit optimism,", m), "table_samefit_optimism_by_method.csv")

# ---------------------------------------------------------------------
# harmonisation and identities
# ---------------------------------------------------------------------
idt <- tab("table_identity_verification.csv")
add("ident_bs_rmse", idt[grepl("1a", result), max_abs_dev], sprintf("%.1e", idt[grepl("1a", result), max_abs_dev]), "max |Brier - RMSE^2|", "table_identity_verification.csv")
add("ident_r2", idt[grepl("1b", result), max_abs_dev], sprintf("%.1e", idt[grepl("1b", result), max_abs_dev]), "max |R2 - (1 - Brier / pi(1 - pi))|", "table_identity_verification.csv")
add("ident_prop1_viol", idt[grepl("Proposition 1", result), n_dev_over_1e3], fc(idt[grepl("Proposition 1", result), n_dev_over_1e3]), "fits violating MAE^2 <= Brier <= MAE", "table_identity_verification.csv")
add("tie_fits", idt[grepl("1c", result), n_fits], fc(idt[grepl("1c", result), n_fits]), "fits with finite Spearman", "table_identity_verification.csv")
add("tie_affected", idt[grepl("1c", result), n_dev_over_1e3], fc(idt[grepl("1c", result), n_dev_over_1e3]), "fits departing from the no tie form by more than 0.001", "table_identity_verification.csv")
add("tie_max", idt[grepl("1c", result), max_abs_dev], f3(idt[grepl("1c", result), max_abs_dev]), "largest departure from the no tie form", "table_identity_verification.csv")
add("tie_median", idt[grepl("1c", result), median_abs_dev], sprintf("%.1e", idt[grepl("1c", result), median_abs_dev]), "median departure from the no tie form", "table_identity_verification.csv")
ds <- dc[is.finite(Spearman)]
ds[, res_1c := { n1 <- round(n_test * test_prev); n0 <- n_test - n1
                 abs(Spearman - (AUC - 0.5) * sqrt(12 * n0 * n1 / (n_test^2 - 1))) }]
ties <- ds[, .(n = .N, n_tie_affected = sum(res_1c > 1e-3), share = mean(res_1c > 1e-3),
               max_dev = signif(max(res_1c), 3)), by = method][order(-share)]
fwrite(ties, file.path(RES, "table_spearman_tie_effect_by_method.csv"))
for (m in METHODS) add(paste0("tie_share_", m), ties[method == m, share], fp(ties[method == m, share], 1), paste("share of fits departing from the no tie form,", m), "table_spearman_tie_effect_by_method.csv")
add("r2_min", min(dc$R2), sprintf("%.0f", min(dc$R2)), "lowest R2 over converged fits", "fit_metrics.rds")
add("r2_max", max(dc$R2), sprintf("%.2f", max(dc$R2)), "highest R2 over converged fits", "fit_metrics.rds")
H <- tab("table_harmonisation_theory.csv")
H[, pair := paste(from, "->", to)]
ratio <- H[!to %in% "R2", median(rel_resid_idr / rel_resid_range)]
add("idr_range_factor", ratio, sprintf("%.1f", ratio), "median ratio of interdecile to range normalised residual, targets other than R2", "table_harmonisation_theory.csv")
for (i in seq_len(nrow(H))) {
  k <- gsub(" -> ", "_", H$pair[i])
  add(paste0("h_rr_", k), H$rel_resid_range[i], f3(H$rel_resid_range[i]), paste("relative residual (range),", H$pair[i]), "table_harmonisation_theory.csv")
  add(paste0("h_ri_", k), H$rel_resid_idr[i], f3(H$rel_resid_idr[i]), paste("relative residual (interdecile),", H$pair[i]), "table_harmonisation_theory.csv")
  add(paste0("h_ts_", k), H$target_shift[i], f3(H$target_shift[i]), paste("target shift,", H$pair[i]), "table_harmonisation_theory.csv")
  add(paste0("h_ss_", k), H$survey_shift[i], f3(H$survey_shift[i]), paste("survey shift,", H$pair[i]), "table_harmonisation_theory.csv")
  add(paste0("h_stab_", k), H$stab_range[i], sprintf("%.2f", H$stab_range[i]), paste("bootstrap stability of the range class,", H$pair[i]), "table_harmonisation_theory.csv")
  add(paste0("h_stabidr_", k), H$stab_idr[i], sprintf("%.2f", H$stab_idr[i]), paste("bootstrap stability of the interdecile class,", H$pair[i]), "table_harmonisation_theory.csv")
}
add("pool_exceeds", sum(H$pooled_exceeds_bound), fc(sum(H$pooled_exceeds_bound)), "translations whose pooled residual exceeds the Corollary 3 bound", "table_harmonisation_theory.csv")
add("within_below_pooled", sum(H$within_context_sd < H$resid_sd), fc(sum(H$within_context_sd < H$resid_sd)), "translations whose within context residual is below the pooled residual", "table_harmonisation_theory.csv")
add("survey_binds_16", sum(H$binding_shift == "survey"), fc(sum(H$binding_shift == "survey")), "translations, of 16, for which the survey shift binds", "table_harmonisation_theory.csv")
orig11 <- c("AUC -> TSS", "TSS -> AUC", "TSS -> Spearman", "AUC -> Spearman", "Spearman -> AUC", "AUC -> RMSE",
            "AUC -> Brier", "AUC -> R2", "R2 -> RMSE", "RMSE -> Brier", "Brier -> RMSE")
add("survey_binds_11", H[pair %in% orig11, sum(binding_shift == "survey")], fc(H[pair %in% orig11, sum(binding_shift == "survey")]), "translations, of the original 11, for which the survey shift binds", "table_harmonisation_theory.csv")
add("stable_range_n", sum(H$stab_range == 1), fc(sum(H$stab_range == 1)), "pairs whose range class held in every bootstrap replicate", "table_harmonisation_theory.csv")
add("stable_idr_n", sum(H$stab_idr == 1), fc(sum(H$stab_idr == 1)), "pairs whose interdecile class held in every bootstrap replicate", "table_harmonisation_theory.csv")
as <- tab("table_directional_asymmetry.csv")
for (i in seq_len(nrow(as))) {
  k <- gsub(" -> ", "_", as$forward[i])
  occ <- as$var_over_range2_from[i] / as$var_over_range2_to[i]
  cr  <- (1 - as$eta2_rev[i]) / (1 - as$eta2_fwd[i])
  add(paste0("asym_obs_", k), as$observed_ratio_sq[i], sprintf("%.2f", as$observed_ratio_sq[i]), paste("observed squared ratio of relative residuals,", as$forward[i]), "table_directional_asymmetry.csv")
  add(paste0("asym_occ_", k), occ, sprintf("%.2f", occ), paste("range occupancy factor,", as$forward[i]), "table_directional_asymmetry.csv")
  add(paste0("asym_cr_", k), cr, sprintf("%.2f", cr), paste("correlation ratio factor,", as$forward[i]), "table_directional_asymmetry.csv")
}
ce <- tab("table_calibration_explains_residual.csv")
add("calib_r2", ce$value[1], fp(ce$value[1], 0), "share of AUC to Brier residual explained by calibration slope and prevalence", "table_calibration_explains_residual.csv")
add("calib_r2_int", ce$value[2], fp(ce$value[2], 0), "adding the calibration intercept", "table_calibration_explains_residual.csv")
add("brt_share_fits", ce$value[3], fp(ce$value[3], 1), "BRT share of fits", "table_calibration_explains_residual.csv")
add("brt_share_rss", ce$value[4], fp(ce$value[4], 1), "BRT share of residual sum of squares", "table_calibration_explains_residual.csv")
rf <- tab("table_translation_to_refinement.csv")
for (t in c("REF", "Brier")) {
  add(paste0("refine_range_", t), rf[target == t, rel_resid_range], f3(rf[target == t, rel_resid_range]), paste("AUC to", t, "relative residual, range"), "table_translation_to_refinement.csv")
  add(paste0("refine_idr_", t), rf[target == t, rel_resid_idr], f3(rf[target == t, rel_resid_idr]), paste("AUC to", t, "relative residual, interdecile"), "table_translation_to_refinement.csv")
}
pc <- tab("table_prevalence_conditioning.csv")
for (i in seq_len(nrow(pc))) {
  k <- paste(pc$from[i], pc$to[i], sep = "_")
  add(paste0("prev_ts_pooled_", k), pc$target_shift_pooled[i], f3(pc$target_shift_pooled[i]), paste("target shift pooled,", k), "table_prevalence_conditioning.csv")
  add(paste0("prev_ts_within_", k), pc$target_shift_within_prev[i], f3(pc$target_shift_within_prev[i]), paste("target shift within prevalence deciles,", k), "table_prevalence_conditioning.csv")
  add(paste0("prev_ss_pooled_", k), pc$survey_shift_pooled[i], f3(pc$survey_shift_pooled[i]), paste("survey shift pooled,", k), "table_prevalence_conditioning.csv")
  add(paste0("prev_ss_within_", k), pc$survey_shift_within_prev[i], f3(pc$survey_shift_within_prev[i]), paste("survey shift within prevalence deciles,", k), "table_prevalence_conditioning.csv")
}

# ---------------------------------------------------------------------
# failure modes
# ---------------------------------------------------------------------
cb <- tab("table_calibration_by_method.csv"); of <- tab("table_overfitting_by_method.csv")
for (m in METHODS) {
  add(paste0("calslope_", m), cb[method == m, calib_slope_med], sprintf("%.2f", cb[method == m, calib_slope_med]), paste("median calibration slope,", m), "table_calibration_by_method.csv")
  add(paste0("overconf_", m), cb[method == m, frac_overconfident], fp(cb[method == m, frac_overconfident], 0), paste("share of fits with slope below 0.8,", m), "table_calibration_by_method.csv")
  add(paste0("gap_", m), of[method == m, gap_AUC], f3(of[method == m, gap_AUC]), paste("train minus test AUC,", m), "table_overfitting_by_method.csv")
}

# ---------------------------------------------------------------------
# heterogeneity across surveys
# ---------------------------------------------------------------------
op <- tab("table_optimism_by_dataset_pooled.csv"); ob <- tab("table_optimism_by_dataset.csv")
for (s in op$dataset) {
  k <- gsub("-", "", s)
  add(paste0("optds_", k), op[dataset == s, optimism], f3(op[dataset == s, optimism]), paste("optimism,", s), "table_optimism_by_dataset_pooled.csv")
  add(paste0("optds_se_", k), op[dataset == s, opt_se], f4(op[dataset == s, opt_se]), paste("its standard error,", s), "table_optimism_by_dataset_pooled.csv")
  ex <- ob[dataset == s & spatial == TRUE, optimism] - ob[dataset == s & spatial == FALSE, optimism]
  add(paste0("excess_", k), ex, f3(ex), paste("spatial minus non spatial optimism,", s), "table_optimism_by_dataset.csv")
  ex2 <- w[dataset == s & method %in% SM, mean(opt)] - w[dataset == s & spatial == FALSE, mean(opt)]
  add(paste0("excess_exJ_", k), ex2, f3(ex2), paste("three spatial smoothers minus non spatial optimism,", s), "fit_metrics.rds")
}
for (s in op$dataset) add(paste0("condauc_", gsub("-", "", s)), op[dataset == s, conditional], f3(op[dataset == s, conditional]),
                          paste("mean conditional AUC over paired draws,", s), "table_optimism_by_dataset_pooled.csv")
dd <- op[dataset == "EBS", optimism] - op[dataset == "NEUS", optimism]
add("optds_EBS_minus_NEUS", dd, f3(dd), "optimism, EBS minus NEUS", "table_optimism_by_dataset_pooled.csv")
pr <- tab("table_survey_profiles.csv")
add("hps_EBS", pr[dataset == "EBS", hauls_per_station], sprintf("%.2f", pr[dataset == "EBS", hauls_per_station]), "hauls per unique station, EBS", "table_survey_profiles.csv")
add("hps_min_other", pr[dataset != "EBS", min(hauls_per_station)], sprintf("%.2f", pr[dataset != "EBS", min(hauls_per_station)]), "hauls per unique station, lowest of the other surveys", "table_survey_profiles.csv")
add("hps_max_other", pr[dataset != "EBS", max(hauls_per_station)], sprintf("%.2f", pr[dataset != "EBS", max(hauls_per_station)]), "hauls per unique station, highest of the other surveys", "table_survey_profiles.csv")
add("optds_ratio", max(op$optimism) / min(op$optimism), sprintf("%.1f", max(op$optimism) / min(op$optimism)), "largest over smallest survey optimism", "table_optimism_by_dataset_pooled.csv")
sd_ <- tab("table_spatial_effect_by_dataset.csv")
for (i in seq_len(nrow(sd_))) {
  k <- paste0(gsub("-", "", sd_$dataset[i]), "_", substr(sd_$target[i], 1, 4))
  add(paste0("spadvds_", k), sd_$spatial_advantage[i], fs(sd_$spatial_advantage[i]), paste("spatial advantage,", sd_$dataset[i], sd_$target[i]), "table_spatial_effect_by_dataset.csv")
  v <- d0[dataset == sd_$dataset[i] & target == sd_$target[i] & method %in% SM, mean(AUC, na.rm = TRUE)] -
       d0[dataset == sd_$dataset[i] & target == sd_$target[i] & spatial == FALSE, mean(AUC, na.rm = TRUE)]
  add(paste0("spadvds_exJ_", k), v, fs(v), paste("spatial advantage without JointSDM,", sd_$dataset[i], sd_$target[i]), "fit_metrics.rds")
}
add("spadvds_negative", sum(sd_$spatial_advantage < 0), fc(sum(sd_$spatial_advantage < 0)), "survey by target cells with a negative spatial advantage", "table_spatial_effect_by_dataset.csv")
cd_ <- tab("table_coverage_by_dataset.csv")
for (s in cd_$dataset) add(paste0("covlift_", gsub("-", "", s)), cd_[dataset == s, coverage_lift], f3(cd_[dataset == s, coverage_lift]), paste("coverage lift,", s), "table_coverage_by_dataset.csv")

# ---------------------------------------------------------------------
# tests of the theoretical predictions
# ---------------------------------------------------------------------
bn <- tab("table_brier_optimism_by_n.csv")
for (sp in c(TRUE, FALSE)) for (nn in SAMPLE_SIZES)
  add(sprintf("bopt_%s_n%d", ifelse(sp, "sp", "ns"), nn), bn[spatial == sp & n == nn, mean_opt_brier], f3(bn[spatial == sp & n == nn, mean_opt_brier]),
      sprintf("Brier optimism, %s methods, n = %d", ifelse(sp, "spatial", "non spatial"), nn), "table_brier_optimism_by_n.csv")
clus_fit <- function(fit, cl, name, src) {
  co <- coef(fit); V <- vcovCL(fit, cluster = cl)
  b_ns <- co["log10n"]; b_int <- co["log10n:spatialTRUE"]
  se_ns <- sqrt(V["log10n", "log10n"]); se_int <- sqrt(V["log10n:spatialTRUE", "log10n:spatialTRUE"])
  b_sp <- b_ns + b_int
  se_sp <- sqrt(V["log10n", "log10n"] + V["log10n:spatialTRUE", "log10n:spatialTRUE"] + 2 * V["log10n", "log10n:spatialTRUE"])
  ols <- summary(fit)$coefficients
  add(paste0(name, "_slope_ns"), b_ns, f4(b_ns), "slope on log10 n, non spatial methods", src)
  add(paste0(name, "_slope_sp"), b_sp, f4(b_sp), "slope on log10 n, spatial methods", src)
  add(paste0(name, "_int"), b_int, f4(b_int), "interaction (spatial minus non spatial slope)", src)
  add(paste0(name, "_int_se_ols"), ols["log10n:spatialTRUE", 2], f4(ols["log10n:spatialTRUE", 2]), "interaction standard error, ordinary least squares", src)
  add(paste0(name, "_slope_ns_se_cl"), se_ns, f4(se_ns), "non spatial slope standard error, clustered by species", src)
  add(paste0(name, "_slope_sp_se_cl"), se_sp, f4(se_sp), "spatial slope standard error, clustered by species", src)
  add(paste0(name, "_int_se_cl"), se_int, f4(se_int), "interaction standard error, clustered by species", src)
  p <- 2 * pnorm(-abs(b_int / se_int))
  add(paste0(name, "_int_p_cl"), p, sprintf("%.1e", p), "interaction p value, clustered by species, normal reference", src)
}
wb <- dcast(dc, dataset + sp_id + n + coverage + rep + method + spatial ~ target, value.var = "Brier")
wb <- wb[is.finite(conditional) & is.finite(marginal)][, `:=`(opt_brier = marginal - conditional, log10n = log10(n))]
clus_fit(lm(opt_brier ~ log10n * spatial, data = wb), ~sp_id, "cor5", "fit_metrics.rds (matched draws, as table_brier_optimism_slope.csv)")
mm[, `:=`(sf_brier = Brier - interp_Brier, log10n = log10(n))]
clus_fit(lm(sf_brier ~ log10n * spatial, data = mm), ~sp_id, "cor5sf", "fit_metrics.rds (same fit, as table_corollary5_samefit.csv)")
wbx <- wb[method != "JointSDM"]
clus_fit(lm(opt_brier ~ log10n * spatial, data = wbx), ~sp_id, "cor5exJ", "fit_metrics.rds (matched draws, JointSDM excluded)")
# the spatial basis dimension is capped below 100 hauls (R/03_methods.R), so the slopes are
# repeated on 100, 200 and 500 hauls, where every smoother has its design basis
clus_fit(lm(opt_brier ~ log10n * spatial, data = wb[n >= 100]), ~sp_id, "cor5n100", "fit_metrics.rds (matched draws, n >= 100)")
clus_fit(lm(sf_brier ~ log10n * spatial, data = mm[n >= 100]), ~sp_id, "cor5sfn100", "fit_metrics.rds (same fit, n >= 100)")
lm5 <- tab("table_lemma5_covariate_count.csv")
add("lemma5_rho", lm5$value[1], sprintf("%.2f", lm5$value[1]), "Spearman, marginal spatial advantage against covariate count", "table_lemma5_covariate_count.csv")
add("lemma5_ns", lm5$value[3], fs(lm5$value[3]), "mean marginal spatial advantage, NS-IBTS", "table_lemma5_covariate_count.csv")
add("lemma5_other", lm5$value[4], fs(lm5$value[4]), "mean marginal spatial advantage, three covariate surveys", "table_lemma5_covariate_count.csv")
n_env <- c(EBS = 3L, `NS-IBTS` = 1L, GMEX = 3L, NEUS = 3L)
ax <- d0[target == "marginal" & method != "JointSDM", .(AUC = mean(AUC, na.rm = TRUE)), by = .(dataset, species, spatial)]
ax <- dcast(ax, dataset + species ~ spatial, value.var = "AUC"); setnames(ax, c("FALSE", "TRUE"), c("ns", "sp"))
ax[, adv := sp - ns][, nenv := n_env[dataset]]
add("lemma5_rho_exJ", cor(ax$adv, ax$nenv, method = "spearman"), sprintf("%.2f", cor(ax$adv, ax$nenv, method = "spearman")), "as lemma5_rho, JointSDM excluded", "fit_metrics.rds")
add("lemma5_ns_exJ", ax[dataset == "NS-IBTS", mean(adv)], fs(ax[dataset == "NS-IBTS", mean(adv)]), "as lemma5_ns, JointSDM excluded", "fit_metrics.rds")
add("lemma5_other_exJ", ax[dataset != "NS-IBTS", mean(adv)], fs(ax[dataset != "NS-IBTS", mean(adv)]), "as lemma5_other, JointSDM excluded", "fit_metrics.rds")
mb <- tab("table_mae_minus_2brier_by_method.csv")
for (m in METHODS) {
  add(paste0("mae2b_", m), mb[method == m, median_mae_minus_2brier], f3(mb[method == m, median_mae_minus_2brier]), paste("median MAE minus 2 Brier,", m), "table_mae_minus_2brier_by_method.csv")
  add(paste0("mae2b_rho_", m), mb[method == m, spearman_with_abs_log_slope], sprintf("%.2f", mb[method == m, spearman_with_abs_log_slope]), paste("Spearman with |log calibration slope|,", m), "table_mae_minus_2brier_by_method.csv")
}
dc[, `:=`(m2b = MAE - 2 * Brier, als = abs(log(pmax(calib_slope, 1e-3))))]
rho_all <- dc[is.finite(calib_slope), cor(m2b, als, method = "spearman")]
add("mae2b_rho_all", rho_all, sprintf("%.2f", rho_all), "Spearman of MAE minus 2 Brier with |log calibration slope|, all fits", "fit_metrics.rds")
od <- tab("table_optimism_decomposition.csv")
for (s in od$dataset) {
  k <- gsub("-", "", s)
  add(paste0("intr_", k), od[dataset == s, intrinsic], fs(od[dataset == s, intrinsic]), paste("intrinsic term (near oracle Brier, block minus extent),", s), "table_optimism_decomposition.csv")
  add(paste0("est_sp_", k), od[dataset == s, estimation_spatial], f3(od[dataset == s, estimation_spatial]), paste("estimation term, spatial methods,", s), "table_optimism_decomposition.csv")
  add(paste0("est_ns_", k), od[dataset == s, estimation_nonspatial], f3(od[dataset == s, estimation_nonspatial]), paste("estimation term, non spatial methods,", s), "table_optimism_decomposition.csv")
  add(paste0("sfbopt_", k), od[dataset == s, opt_brier], f3(od[dataset == s, opt_brier]), paste("same fit Brier optimism,", s), "table_optimism_decomposition.csv")
}
add("intr_var_ratio", var(od$intrinsic) / var(od$opt_brier), sprintf("%.2f", var(od$intrinsic) / var(od$opt_brier)), "variance of intrinsic term over variance of Brier optimism, across surveys", "table_optimism_decomposition.csv")
c4 <- tab("table_corollary4_regression.csv")
for (t in c4$term) {
  k <- gsub("[^A-Za-z0-9]", "", t)
  add(paste0("c4_", k), c4[term == t, estimate], f3(c4[term == t, estimate]), paste("Corollary 4 coefficient,", t), "table_corollary4_regression.csv")
  add(paste0("c4se_", k), c4[term == t, se], f4(c4[term == t, se]), paste("its standard error,", t), "table_corollary4_regression.csv")
}
add("c4_nspecies", c4$n_species[1], fc(c4$n_species[1]), "species in the Corollary 4 model (identified residual range)", "table_corollary4_regression.csv")
add("c4_nfits", c4$n_fits[1], fc(c4$n_fits[1]), "fits in the Corollary 4 model (identified residual range)", "table_corollary4_regression.csv")
c4s <- tab("table_corollary4_sensitivity.csv")
for (sp in unique(c4s$specification)) {
  pre <- if (startsWith(sp, "original")) "c4orig_" else "c4km_"
  for (t in c4s[specification == sp, term]) {
    k <- gsub("[^A-Za-z0-9]", "", t)
    add(paste0(pre, k), c4s[specification == sp & term == t, estimate], f3(c4s[specification == sp & term == t, estimate]), paste0("Corollary 4 coefficient, ", sp, ", ", t), "table_corollary4_sensitivity.csv")
    add(paste0(pre, "se_", k), c4s[specification == sp & term == t, se], f4(c4s[specification == sp & term == t, se]), paste0("its standard error, ", sp, ", ", t), "table_corollary4_sensitivity.csv")
  }
}
# residual correlation range regimes (05d_practical_range.R) and block separation (06b, part G)
pr <- tab("practical_range.csv")
for (g in c("identified", "no structure", "below lags", "no sill"))
  add(paste0("range_regime_", gsub(" ", "", g)), pr[regime == g, .N], fc(pr[regime == g, .N]), paste("species whose residual variogram regime is", g), "practical_range.csv")
add("range_not_identified", pr[regime != "identified", .N], fc(pr[regime != "identified", .N]), "species whose residual variogram does not identify a range", "practical_range.csv")
add("range_identified_median", pr[regime == "identified", median(range_km)], fc(pr[regime == "identified", median(range_km)]), "median identified practical range (km)", "practical_range.csv")
add("range_identified_min", pr[regime == "identified", min(range_km)], fc(pr[regime == "identified", min(range_km)]), "smallest identified practical range (km)", "practical_range.csv")
add("range_identified_max", pr[regime == "identified", max(range_km)], fc(pr[regime == "identified", max(range_km)]), "largest identified practical range (km)", "practical_range.csv")
bs_ <- tab("table_block_separation.csv")[dataset == "all"]
add("sep_draws", bs_$draws, fc(bs_$draws), "marginal draws with a recorded block separation", "table_block_separation.csv")
add("sep_median_km", bs_$median_separation_km, fc(bs_$median_separation_km), "median separation of the block from its training locations (km)", "table_block_separation.csv")
add("sep_draws_identified", bs_$draws_identified, fc(bs_$draws_identified), "marginal draws of species with an identified residual range", "table_block_separation.csv")
add("sep_scaled_median_identified", bs_$median_scaled_identified, sprintf("%.2f", bs_$median_scaled_identified), "median separation in identified practical ranges", "table_block_separation.csv")
add("sep_share_beyond_identified", bs_$share_beyond_one_range_identified, fp(bs_$share_beyond_one_range_identified, 0), "share of those draws with the block one practical range or more away", "table_block_separation.csv")
add("matern32_corr_at_range", 2 * exp(-1), sprintf("%.2f", 2 * exp(-1)), "Matern (nu = 3/2) correlation at a distance equal to its range parameter, (1 + 1) exp(-1)", "analytic")
fd <- tab("table_fill_distance_model.csv")
add("fill_coef", fd[term == "log2_fill", estimate], sprintf("%+.4f", fd[term == "log2_fill", estimate]), "marginal AUC per doubling of fill distance", "table_fill_distance_model.csv")
add("fill_se", fd[term == "log2_fill", se], f4(fd[term == "log2_fill", se]), "its standard error", "table_fill_distance_model.csv")
add("corollary2_limit", nrow(dc)^-0.5, sprintf("%.4f", nrow(dc)^-0.5), "N^(-1/2) for the converged fit count", "fit_metrics.rds")

out <- rbindlist(N); setnames(out, "id", "key")
fwrite(out, file.path(RES, "manuscript_numbers.csv"))
cat(sprintf("wrote %d numbers to results/manuscript_numbers.csv\n", nrow(out)))
