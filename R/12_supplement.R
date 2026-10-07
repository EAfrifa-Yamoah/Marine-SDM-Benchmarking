# =====================================================================
# 12_supplement.R
# Tables and figures of the Supporting Information, from the full run.
#   results/supplement/table_S*.csv
#   figures/supplement/fig_S*.png
# Colour follows R/palette.R (hue identifies a method and nothing else).
# Run after 06 to 11. Object names inside the script (S6, S7, ...) follow the order of
# computation; the output file names carry the numbering used in the Supporting Information.
# =====================================================================
suppressMessages({ library(data.table); library(ggplot2); library(patchwork); library(viridis) })
if (!nzchar(Sys.getenv("SMD_ROOT"))) Sys.setenv(SMD_ROOT = getwd())
setwd(Sys.getenv("SMD_ROOT")); source(file.path("R", "config.R")); source(file.path("R", "palette.R"))
RES <- PATHS$results; FIG <- PATHS$figures
ST <- file.path(RES, "supplement"); SF <- file.path(FIG, "supplement")
dir.create(ST, showWarnings = FALSE, recursive = TRUE); dir.create(SF, showWarnings = FALSE, recursive = TRUE)
tab <- function(f) fread(file.path(RES, f))
wt <- function(x, name) { fwrite(x, file.path(ST, name)); cat("wrote", name, "\n") }
r3 <- function(x) round(x, 3)
SURVEY <- c(EBS = "Eastern Bering Sea", `NS-IBTS` = "North Sea (IBTS)", GMEX = "Gulf of Mexico", NEUS = "Northeast US shelf")
SM <- c("SpatialRF", "SpatialGAM", "GeostatGP")
theme_s <- theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), plot.title = element_blank(),
        strip.text = element_text(face = "bold"))
savef <- function(p, name, w = 8, h = 5) { ggsave(file.path(SF, name), p, width = w, height = h, dpi = 200, bg = "white"); cat("wrote", name, "\n") }

d <- as.data.table(readRDS(file.path(RES, "fit_metrics.rds")))
d[, method := factor(method, levels = METHODS)]
d[, sp_id := paste(dataset, species)]

# ---- Table S1 species ---------------------------------------------------
S1 <- rbindlist(lapply(names(SURVEY), function(code) {
  r <- readRDS(file.path(PATHS$processed, paste0(code, "_processed.rds")))
  s <- as.data.table(r$species)
  s[, .(survey = SURVEY[code], species = accepted_name, band, presences = n_pres,
        prevalence_pct = round(100 * prevalence, 1), niche_breadth = round(niche, 2), trophic_proxy = trophic)]
}))
S1[, band := factor(band, levels = c("rare", "intermediate", "common"))]
setorder(S1, survey, band, -prevalence_pct)
wt(S1, "table_S1_species.csv")

# ---- Table S2 survey profiles -------------------------------------------
pr <- tab("table_survey_profiles.csv")
S2 <- pr[, .(survey = SURVEY[dataset], hauls, years, unique_stations, hauls_per_station,
             extent_km2, span_x_km, span_y_km, covariates = gsub("\\+", ", ", env))]
wt(S2, "table_S2_survey_profiles.csv")

# ---- Table S3 degenerate splits -----------------------------------------
sf <- tab("split_failures.csv")
attempted <- 20 * 3 * 2 * N_REPLICATES       # species per band x coverages x targets x replicates, per n
S3 <- sf[, .N, by = .(band, n)]
S3 <- dcast(S3, band ~ n, value.var = "N", fill = 0)
S3 <- rbind(S3, data.table(band = "common", `30` = 0, `50` = 0, `100` = 0, `200` = 0, `500` = 0), fill = TRUE)
S3[, total := rowSums(.SD), .SDcols = as.character(SAMPLE_SIZES)]
S3[, rate_at_n30_pct := round(100 * `30` / attempted, 1)]
S3[, band := factor(band, levels = c("rare", "intermediate", "common"))]; setorder(S3, band)
wt(S3, "table_S3_degenerate_splits.csv")

# ---- Table S4 performance by method and target ----------------------------
S4 <- d[, .(AUC = r3(mean(AUC, na.rm = TRUE)), TSS = r3(mean(TSS, na.rm = TRUE)),
            Spearman = r3(mean(Spearman, na.rm = TRUE)), Brier = r3(mean(Brier, na.rm = TRUE)),
            RMSE = r3(mean(RMSE, na.rm = TRUE)), R2_median = r3(median(R2, na.rm = TRUE)),
            calibration_slope_median = round(median(calib_slope, na.rm = TRUE), 2),
            train_minus_test_AUC = r3(mean(gap_AUC, na.rm = TRUE)), fits = .N), by = .(target, method)]
setorder(S4, target, method)
wt(S4, "table_S4_performance_by_method_target.csv")

# ---- Table S5 AUC by method, target and sample size -----------------------
a5 <- d[, .(AUC = mean(AUC, na.rm = TRUE)), by = .(target, method, n)]
S5 <- dcast(a5, method ~ target + n, value.var = "AUC")
for (cc in setdiff(names(S5), "method")) set(S5, j = cc, value = r3(S5[[cc]]))
sl <- tab("table_scaling_slopes.csv")[, .(method, gain_per_decade_conditional = r3(slope_per_decade))]
S5 <- merge(S5, sl, by = "method"); S5[, method := factor(method, levels = METHODS)]; setorder(S5, method)
wt(S5, "table_S5_auc_by_method_target_n.csv")

# ---- Table S7 harmonisation, both rules -----------------------------------
H <- tab("table_harmonisation_theory.csv")
S6 <- H[, .(translation = paste(from, "→", to), kind, pooled_residual_sd = signif(resid_sd, 3),
            within_context_sd = signif(within_context_sd, 3), corollary3_bound = signif(pool_bound, 3),
            binding_shift, class_range, stability_range = stab_range, class_interdecile = class_idr,
            stability_interdecile = stab_idr)]
wt(S6, "table_S7_harmonisation_both_rules.csv")

# ---- Table S6 hierarchical model ------------------------------------------
fe <- tab("table_fixed_effects.csv"); vp <- tab("table_variance_partition.csv")
S7 <- rbind(fe[, .(component = "fixed effect", term, estimate = signif(Estimate, 3), se = signif(`Std. Error`, 2), share_of_variance = NA_real_)],
            vp[, .(component = "variance component", term = grp, estimate = signif(vcov, 3), se = NA_real_, share_of_variance = round(prop, 3))])
wt(S7, "table_S6_hierarchical_model.csv")

# ---- Table S9 optimism by method -----------------------------------------
w <- dcast(d, dataset + sp_id + band + n + coverage + rep + method + spatial ~ target, value.var = "AUC")
w <- w[!is.na(conditional) & !is.na(marginal)][, opt := conditional - marginal]
spm <- w[, .(o = mean(opt)), by = .(method, sp_id)][, .(se_species = sd(o) / sqrt(.N)), by = method]
sfm <- tab("table_samefit_optimism_by_method.csv")
vo <- tab("table_validation_optimism.csv")
S8 <- merge(vo[, .(method, spatial, conditional_AUC = r3(cond), marginal_AUC = r3(marg), optimism = r3(optimism),
                   se_naive = signif(opt_sd / sqrt(n_pairs), 2), pairs = n_pairs)],
            spm[, .(method = as.character(method), se_species = signif(se_species, 2))], by = "method")
S8 <- merge(S8, sfm[, .(method, samefit_optimism_AUC = r3(opt_AUC), samefit_optimism_Brier = r3(opt_Brier))], by = "method")
setorder(S8, -optimism)
wt(S8, "table_S9_optimism_by_method.csv")

# ---- Table S10 optimism by sample size, coverage and class --------------------
w[, cls := ifelse(method %in% SM, "spatial smoothers", ifelse(method == "JointSDM", "joint model", "non-spatial"))]
S9 <- dcast(w[, .(o = r3(mean(opt))), by = .(n, cls)], n ~ cls, value.var = "o")
S9b <- dcast(w[, .(o = r3(mean(opt))), by = .(coverage, cls)], coverage ~ cls, value.var = "o")
setnames(S9, "n", "level"); setnames(S9b, "coverage", "level")
S9 <- rbind(S9[, factor := "training sample size"][, level := as.character(level)], S9b[, factor := "spatial coverage"])
setcolorder(S9, c("factor", "level"))
wt(S9, "table_S10_optimism_by_n_and_coverage.csv")

# ---- Table S11 survey heterogeneity ---------------------------------------
ob <- tab("table_optimism_by_dataset.csv"); op <- tab("table_optimism_by_dataset_pooled.csv")
sps <- w[, .(o = mean(opt)), by = .(dataset, sp_id)][, .(se_species = sd(o) / sqrt(.N)), by = dataset]
sa <- d[, .(AUC = mean(AUC, na.rm = TRUE)), by = .(dataset, target, cls = ifelse(method %in% SM, "sm", ifelse(spatial, "jsdm", "ns")))]
adv <- function(ds, tg, with_j) {
  x <- sa[dataset == ds & target == tg]
  sp <- if (with_j) d[dataset == ds & target == tg & spatial == TRUE, mean(AUC, na.rm = TRUE)] else x[cls == "sm", AUC]
  sp - x[cls == "ns", AUC]
}
S10 <- op[, .(survey = SURVEY[dataset], dataset, optimism = r3(optimism))]
S10 <- merge(S10, sps[, .(dataset, se_species = signif(se_species, 2))], by = "dataset")
S10[, `:=`(excess_spatial = r3(ob[dataset == .BY$dataset & spatial == TRUE, optimism] - ob[dataset == .BY$dataset & spatial == FALSE, optimism])), by = dataset]
S10[, `:=`(excess_spatial_without_joint = r3(w[dataset == .BY$dataset & method %in% SM, mean(opt)] - w[dataset == .BY$dataset & spatial == FALSE, mean(opt)]),
           advantage_conditional = r3(adv(.BY$dataset, "conditional", TRUE)), advantage_conditional_without_joint = r3(adv(.BY$dataset, "conditional", FALSE)),
           advantage_marginal = r3(adv(.BY$dataset, "marginal", TRUE)), advantage_marginal_without_joint = r3(adv(.BY$dataset, "marginal", FALSE))), by = dataset]
S10 <- merge(S10, tab("table_coverage_by_dataset.csv")[, .(dataset, coverage_lift = r3(coverage_lift))], by = "dataset")
setorder(S10, -optimism); S10[, dataset := NULL]
wt(S10, "table_S11_survey_heterogeneity.csv")

# ---- Tables S12 to S14 theory tests: Corollary 5, Corollary 4, fill distance, decomposition
N <- tab("manuscript_numbers.csv"); g <- function(k) N[key == k, value]

S11 <- rbindlist(lapply(names(c(cor5 = 1, cor5sf = 1, cor5exJ = 1)), function(k) data.table(
  contrast = c(cor5 = "matched draws", cor5sf = "same fit", cor5exJ = "matched draws, joint model excluded")[k],
  slope_nonspatial = signif(g(paste0(k, "_slope_ns")), 3), se_nonspatial_clustered = signif(g(paste0(k, "_slope_ns_se_cl")), 2),
  slope_spatial = signif(g(paste0(k, "_slope_sp")), 3), se_spatial_clustered = signif(g(paste0(k, "_slope_sp_se_cl")), 2),
  difference = signif(g(paste0(k, "_int")), 3), se_difference_ols = signif(g(paste0(k, "_int_se_ols")), 2),
  se_difference_clustered = signif(g(paste0(k, "_int_se_cl")), 2), p_clustered = signif(g(paste0(k, "_int_p_cl")), 2))))
wt(S11, "table_S12_corollary5_slopes.csv")
c4 <- tab("table_corollary4_regression.csv"); fd <- tab("table_fill_distance_model.csv"); od <- tab("table_optimism_decomposition.csv")
S12 <- rbind(c4[, .(model = "Corollary 4: estimation term", term, estimate = signif(estimate, 3), se = signif(se, 2))],
             fd[, .(model = "Proposition 5: marginal AUC on fill distance", term, estimate = signif(estimate, 3), se = signif(se, 2))])
wt(S12, "table_S13_corollary4_and_fill_distance.csv")
S13 <- od[, .(survey = SURVEY[dataset], intrinsic = r3(intrinsic), estimation_spatial = r3(estimation_spatial),
              estimation_nonspatial = r3(estimation_nonspatial), samefit_brier_optimism = r3(opt_brier))]
wt(S13, "table_S14_optimism_decomposition.csv")

# ---- Table S15 diagnostics ----------------------------------------------
idt <- tab("table_identity_verification.csv"); ties <- tab("table_spearman_tie_effect_by_method.csv")
mb <- tab("table_mae_minus_2brier_by_method.csv"); pc <- tab("table_prevalence_conditioning.csv")
as <- tab("table_directional_asymmetry.csv"); rf <- tab("table_translation_to_refinement.csv")
S14a <- idt[, .(check = result, fits = n_fits, max_abs_dev = signif(max_abs_dev, 2), median_abs_dev = signif(median_abs_dev, 2), fits_over_0.001 = n_dev_over_1e3)]
S14b <- merge(ties[, .(method, share_tie_affected = round(share, 3), max_dev)],
              mb[, .(method, median_mae_minus_2brier = r3(median_mae_minus_2brier), median_calibration_slope = round(median_calib_slope, 2),
                     spearman_with_abs_log_slope = round(spearman_with_abs_log_slope, 2))], by = "method")
S14c <- pc[, .(translation = paste(from, "→", to), target_shift_pooled = r3(target_shift_pooled), target_shift_within_deciles = r3(target_shift_within_prev),
               survey_shift_pooled = r3(survey_shift_pooled), survey_shift_within_deciles = r3(survey_shift_within_prev))]
S14d <- as[, .(pair = paste(forward, "/", reverse), observed_ratio_sq = round(observed_ratio_sq, 2), predicted_ratio_sq = round(predicted_ratio_sq, 2),
               range_occupancy_factor = round(var_over_range2_from / var_over_range2_to, 2), correlation_ratio_factor = round((1 - eta2_rev) / (1 - eta2_fwd), 2))]
wt(S14a, "table_S15a_identities.csv"); wt(S14b, "table_S15b_ties_and_calibration_diagnostic.csv")
wt(S14c, "table_S15c_prevalence_conditioning.csv"); wt(S14d, "table_S15d_directional_asymmetry.csv")
wt(rf[, .(target, rel_resid_range = r3(rel_resid_range), rel_resid_interdecile = r3(rel_resid_idr))], "table_S15e_refinement_translation.csv")

# ---- Table S8 bootstrap class frequencies -------------------------------
bc <- rbindlist(readRDS(file.path(RES, "boot_classes.rds")), idcol = "rep")
S15 <- bc[, .(replicates = .N,
              range_reliable = sum(class_range == "reliable"), range_context = sum(class_range == "context dependent"), range_unreliable = sum(class_range == "unreliable"),
              idr_reliable = sum(class_idr == "reliable"), idr_context = sum(class_idr == "context dependent"), idr_unreliable = sum(class_idr == "unreliable")),
          by = .(translation = paste(from, "→", to))]
wt(S15, "table_S8_bootstrap_class_frequencies.csv")

# ---- Table S16 pilot against full study ----------------------------------
P <- function(f) fread(file.path("results", "pilot", f))
pd <- P("fit_metrics.csv")
pvo <- P("table_validation_optimism.csv"); pse <- P("table_spatial_effect.csv"); pvp <- P("table_variance_partition.csv")
pcov <- pd[target == "marginal", .(a = mean(AUC)), by = coverage]
fcov <- d[target == "marginal", .(a = mean(AUC, na.rm = TRUE)), by = coverage]
pw <- dcast(pd, dataset + species + n + coverage + rep + method + spatial ~ target, value.var = "AUC")[!is.na(conditional) & !is.na(marginal)][, o := conditional - marginal]
pop <- pw[, .(o = mean(o)), by = dataset]
row <- function(q, p, f) data.table(quantity = q, pilot = r3(p), full = r3(f))
S16 <- rbind(
  row("conditional AUC, n = 30 (all methods)", pd[target == "conditional" & n == 30, mean(AUC)], d[target == "conditional" & n == 30, mean(AUC, na.rm = TRUE)]),
  row("conditional AUC, n = 500 (all methods)", pd[target == "conditional" & n == 500, mean(AUC)], d[target == "conditional" & n == 500, mean(AUC, na.rm = TRUE)]),
  row("JointSDM conditional AUC", pd[target == "conditional" & method == "JointSDM", mean(AUC)], d[target == "conditional" & method == "JointSDM", mean(AUC, na.rm = TRUE)]),
  row("spatial advantage, conditional (as specified)", pse[target == "conditional", spatial_advantage], g("spadv_cond")),
  row("spatial advantage, marginal (as specified)", pse[target == "marginal", spatial_advantage], g("spadv_marg")),
  row("optimism, spatial class (as specified)", pw[spatial == TRUE, mean(o)], g("opt_spatial")),
  row("optimism, non-spatial class", pw[spatial == FALSE, mean(o)], g("opt_nonspatial")),
  row("optimism, Eastern Bering Sea", pop[dataset == "EBS", o], g("optds_EBS")),
  row("optimism, North Sea", pop[dataset == "NS-IBTS", o], g("optds_NSIBTS")),
  row("optimism, Gulf of Mexico", pop[dataset == "GMEX", o], g("optds_GMEX")),
  row("optimism, Northeast US shelf", pop[dataset == "NEUS", o], g("optds_NEUS")),
  row("marginal AUC, clustered coverage", pcov[coverage == "clustered", a], fcov[coverage == "clustered", a]),
  row("marginal AUC, distributed coverage", pcov[coverage == "distributed", a], fcov[coverage == "distributed", a]),
  row("variance share, residual", pvp[grp == "Residual", prop], g("vp_residual")),
  row("variance share, species", pvp[grp == "sp_id", prop], g("vp_species")),
  row("variance share, survey", pvp[grp == "dataset", prop], g("vp_survey")),
  row("variance share, survey by target", pvp[grp == "dataset:target", prop], g("vp_survey_target")),
  row("Lemma 5 correlation", P("table_lemma5_covariate_count.csv")$value[1], g("lemma5_rho")))
wt(S16, "table_S16_pilot_against_full.csv")

# ======================================================================
# figures
# ======================================================================
mcol <- METHOD_COLOURS[METHODS]
# S1 AUC against n by survey and target
f1 <- d[, .(AUC = mean(AUC, na.rm = TRUE)), by = .(dataset, target, method, n)]
f1[, dataset := factor(SURVEY[dataset], levels = SURVEY)]
f1[, target := factor(target, c("conditional", "marginal"), c("conditional (interpolation)", "marginal (extrapolation)"))]
p1 <- ggplot(f1, aes(n, AUC, colour = method)) + geom_line(linewidth = 0.7) + geom_point(size = 1.1) +
  facet_grid(target ~ dataset) + scale_x_log10(breaks = SAMPLE_SIZES) + scale_colour_manual(values = mcol) +
  labs(x = "training sample size (log scale)", y = "mean test AUC", colour = "method") + theme_s +
  theme(legend.position = "bottom")
savef(p1, "fig_S1_auc_by_survey.png", w = 11, h = 6.2)

# S3 optimism against n: AUC by method, Brier by class
f2a <- w[, .(o = mean(opt)), by = .(method, n)]
p2a <- ggplot(f2a, aes(n, o, colour = method)) + geom_line(linewidth = 0.7) + geom_point(size = 1.2) +
  scale_x_log10(breaks = SAMPLE_SIZES) + scale_colour_manual(values = mcol) +
  labs(x = "training sample size (log scale)", y = "optimism (AUC)", colour = "method") + theme_s
bn <- tab("table_brier_optimism_by_n.csv")[, cls := factor(ifelse(spatial, "spatial method", "non-spatial method"), levels = names(CLASS_SHAPES))]
p2b <- ggplot(bn, aes(n, mean_opt_brier, shape = cls, linetype = cls)) + geom_line(colour = INK, linewidth = 0.5) +
  geom_point(colour = INK, fill = "white", size = 2.4, stroke = 0.9) + scale_x_log10(breaks = SAMPLE_SIZES) +
  scale_shape_manual(values = CLASS_SHAPES, name = NULL) + scale_linetype_manual(values = c("spatial method" = "solid", "non-spatial method" = "22"), name = NULL) +
  labs(x = "training sample size (log scale)", y = "optimism (Brier score)") + theme_s + theme(legend.position = "bottom")
savef(p2a + p2b + plot_annotation(tag_levels = "a"), "fig_S3_optimism_by_n.png", w = 11, h = 4.6)

# S4 paired against same fit optimism by method
f3 <- merge(S8[, .(method, paired = optimism)], S8[, .(method, samefit = samefit_optimism_AUC)], by = "method")
f3 <- melt(f3, id.vars = "method", variable.name = "contrast", value.name = "o")
f3[, method := factor(method, levels = S8[order(optimism), method])]
f3[, contrast := factor(contrast, c("paired", "samefit"), c("paired draws", "same fit"))]
p3 <- ggplot(f3, aes(o, method, colour = method, shape = contrast)) + geom_line(aes(group = method), colour = "grey70") +
  geom_point(size = 2.8, fill = "white", stroke = 1) + scale_colour_manual(values = mcol, guide = "none") +
  scale_shape_manual(values = c("paired draws" = 16, "same fit" = 21), name = NULL) +
  labs(x = "optimism (conditional or interpolation minus marginal AUC)", y = NULL) + theme_s + theme(legend.position = "bottom")
savef(p3, "fig_S4_paired_against_same_fit.png", w = 7.5, h = 4.6)

# S5 marginal spatial advantage by species, with and without the joint model
adv_sp <- function(with_j) {
  x <- d[target == "marginal" & (with_j | method != "JointSDM"), .(AUC = mean(AUC, na.rm = TRUE)), by = .(dataset, species, spatial)]
  x <- dcast(x, dataset + species ~ spatial, value.var = "AUC"); setnames(x, c("FALSE", "TRUE"), c("ns", "sp"))
  x[, .(dataset, species, advantage = sp - ns, analysis = if (with_j) "as specified" else "joint model excluded")]
}
f4 <- rbind(adv_sp(TRUE), adv_sp(FALSE))
f4[, dataset := factor(SURVEY[dataset], levels = rev(SURVEY))]
p4 <- ggplot(f4, aes(advantage, dataset)) + geom_vline(xintercept = 0, linetype = 2, colour = INK_MID) +
  geom_point(colour = INK, alpha = 0.7, size = 1.8, position = position_jitter(height = 0.12, width = 0, seed = 1)) +
  stat_summary(fun = mean, geom = "point", shape = 124, size = 7, colour = INK) +
  facet_wrap(~ analysis) + labs(x = "marginal spatial advantage per species (spatial minus non-spatial AUC)", y = NULL) + theme_s
savef(p4, "fig_S5_spatial_advantage_by_species.png", w = 9, h = 4.2)

# S2 bootstrap stability of each translation class
f5 <- melt(H[, .(translation = paste(from, "→", to), range = stab_range, interdecile = stab_idr)], id.vars = "translation",
           variable.name = "rule", value.name = "stability")
f5[, translation := factor(translation, levels = H[order(stab_range, stab_idr, decreasing = TRUE), paste(from, "→", to)])]
p5 <- ggplot(f5, aes(stability, translation, shape = rule)) + geom_vline(xintercept = 0.9, linetype = 2, colour = INK_MID) +
  geom_point(colour = INK, size = 2.6, fill = "white", stroke = 0.9, position = position_dodge(width = 0.5)) +
  scale_shape_manual(values = c(range = 16, interdecile = 21), name = "class rule") +
  scale_x_continuous(limits = c(0.8, 1.005), breaks = seq(0.8, 1, 0.05)) +
  labs(x = "share of 100 species bootstrap replicates keeping the full data class", y = NULL) + theme_s
savef(p5, "fig_S2_bootstrap_stability.png", w = 7.5, h = 5.2)

# S6 fill distance by design
fdd <- tab("table_fill_distance_by_design.csv")
fdd[, dataset := factor(SURVEY[dataset], levels = SURVEY)]
fdd[, coverage := factor(coverage, c("clustered", "intermediate", "distributed"))]
p6 <- ggplot(fdd, aes(n, fill_km_mean, linetype = coverage, shape = coverage)) + geom_line(colour = INK) +
  geom_point(colour = INK, fill = "white", size = 2) + facet_wrap(~ dataset, nrow = 1) + scale_x_log10(breaks = SAMPLE_SIZES) +
  scale_shape_manual(values = c(clustered = 16, intermediate = 17, distributed = 21)) +
  labs(x = "training sample size (log scale)", y = "fill distance (km)", linetype = "coverage", shape = "coverage") +
  theme_s + theme(legend.position = "bottom")
savef(p6, "fig_S6_fill_distance.png", w = 11, h = 4)

# S7 estimation term against scaled separation (Corollary 4)
o <- tab("oracle_intrinsic.csv")
ob <- o[, .(bs = median(brier_in_sample), range_km = practical_range_km[1]), by = .(dataset, species, edge, region)]
ow <- dcast(ob, dataset + species + edge + range_km ~ region, value.var = "bs")[, intrinsic := block - extent]
m <- d[converged == TRUE & target == "marginal" & is.finite(interp_Brier) & is.finite(Brier)]
m <- merge(m, ow[, .(dataset, species, edge, intrinsic, range_km)], by = c("dataset", "species", "edge"))
m[, `:=`(estimation = (Brier - interp_Brier) - intrinsic, sep_scaled = separation_km / range_km)]
m <- m[is.finite(estimation) & is.finite(sep_scaled)]
m[, bin := cut(log1p(sep_scaled), breaks = quantile(log1p(sep_scaled), seq(0, 1, 0.1)), include.lowest = TRUE, labels = FALSE)]
f7 <- m[, .(x = median(log1p(sep_scaled)), y = mean(estimation)), by = .(bin, cls = factor(ifelse(spatial, "spatial method", "non-spatial method"), levels = names(CLASS_SHAPES)))]
p7 <- ggplot(f7, aes(x, y, shape = cls, linetype = cls)) + geom_line(colour = INK) + geom_point(colour = INK, fill = "white", size = 2.4, stroke = 0.9) +
  scale_shape_manual(values = CLASS_SHAPES, name = NULL) + scale_linetype_manual(values = c("spatial method" = "solid", "non-spatial method" = "22"), name = NULL) +
  labs(x = "log(1 + separation / practical range), decile medians", y = "estimation term of Brier optimism") + theme_s + theme(legend.position = "bottom")
savef(p7, "fig_S7_estimation_against_separation.png", w = 7, h = 4.6)

# S9 MAE minus 2 Brier against calibration slope
dd <- d[converged == TRUE & is.finite(calib_slope) & calib_slope > 0]
dd[, ls := log10(calib_slope)][, bin := cut(ls, breaks = seq(-2, 1, 0.2), include.lowest = TRUE)]
f8 <- dd[!is.na(bin), .(x = median(calib_slope), y = median(MAE - 2 * Brier), k = .N), by = .(method, bin)]
# bins holding fewer than 500 fits (0.3% of a method's fits) are not drawn; the one
# omitted bin with an extreme value (GeostatGP, slopes 6.3 to 10, 208 fits, median -0.49)
# is reported in the caption of Figure S9
wt(f8[k < 500], "fig_S9_omitted_bins.csv")
f8 <- f8[k >= 500]
p8 <- ggplot(f8, aes(x, y, colour = method)) + geom_hline(yintercept = 0, linetype = 2, colour = INK_MID) + geom_vline(xintercept = 1, linetype = 3, colour = INK_MID) +
  geom_line(linewidth = 0.7) + geom_point(size = 1.4) + scale_x_log10() + scale_colour_manual(values = mcol) +
  labs(x = "calibration slope (log scale), binned", y = "median of MAE minus twice the Brier score", colour = "method") + theme_s
savef(p8, "fig_S9_calibration_diagnostic.png", w = 8, h = 4.8)

# S8 optimism decomposition by survey
f9 <- melt(S13, id.vars = "survey", measure.vars = c("intrinsic", "estimation_spatial", "estimation_nonspatial"), variable.name = "term", value.name = "v")
f9[, term := factor(term, c("intrinsic", "estimation_spatial", "estimation_nonspatial"), c("intrinsic term", "estimation term, spatial methods", "estimation term, non-spatial methods"))]
f9[, survey := factor(survey, levels = rev(SURVEY))]
p9 <- ggplot(f9, aes(v, survey, shape = term)) + geom_vline(xintercept = 0, linetype = 2, colour = INK_MID) +
  geom_point(colour = INK, fill = "white", size = 2.8, stroke = 0.9, position = position_dodge(width = 0.5)) +
  scale_shape_manual(values = c("intrinsic term" = 4, "estimation term, spatial methods" = 16, "estimation term, non-spatial methods" = 21), name = NULL) +
  labs(x = "contribution to same fit Brier optimism", y = NULL) + theme_s + theme(legend.position = "bottom", legend.direction = "vertical")
savef(p9, "fig_S8_optimism_decomposition.png", w = 7, h = 4.8)
cat("supplement complete\n")
