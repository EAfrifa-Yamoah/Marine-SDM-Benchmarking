# =====================================================================
# 08_figures.R
# Publication figures for the benchmark harmonisation pilot.
# Seven figures, one per analytic theme, written as PNG to figures/.
# A consistent minimal theme and the viridis palette are used throughout.
# =====================================================================

suppressMessages({
  library(data.table); library(ggplot2); library(viridis); library(patchwork)
})
source(file.path("R", "config.R"))

RES <- PATHS$results; FIG <- PATHS$figures
d <- as.data.table(readRDS(file.path(RES, "fit_metrics.rds")))
A <- readRDS(file.path(RES, "analysis.rds"))
H <- readRDS(file.path(RES, "harmonise_tables.rds"))

d[, method := factor(method, levels = METHODS)]
d[, target := factor(target, levels = c("conditional", "marginal"),
                     labels = c("conditional (interpolation)",
                                "marginal (extrapolation)"))]
d[, coverage := factor(coverage,
      levels = c("clustered", "intermediate", "distributed"))]
d[, band := factor(band, levels = c("rare", "intermediate", "common"))]

theme_pub <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_blank(),
        plot.subtitle = element_blank(),
        legend.position = "right",
        strip.text = element_text(face = "bold"))
# Okabe-Ito qualitative palette: colour-blind safe across deuteranopia,
# protanopia and tritanopia (orange, sky blue, bluish green, blue,
# vermillion, reddish purple, black)
mcol <- c("#E69F00", "#56B4E9", "#009E73", "#0072B2",
          "#D55E00", "#CC79A7", "#000000")[seq_along(METHODS)]
names(mcol) <- METHODS
savef <- function(p, name, w = 8, h = 5.2) {
  ggsave(file.path(FIG, name), p, width = w, height = h, dpi = 200,
         bg = "white")
  cat("wrote", name, "\n")
}

# ---- Figure 1  sample size scaling (RQ1) ----------------------------
sc <- d[grepl("conditional", target),
        .(AUC = mean(AUC), se = sd(AUC) / sqrt(.N)), by = .(method, n)]
f1 <- ggplot(sc, aes(n, AUC, colour = method)) +
  geom_ribbon(aes(ymin = AUC - se, ymax = AUC + se, fill = method),
              alpha = 0.12, colour = NA) +
  geom_line(linewidth = 0.9) + geom_point(size = 1.6) +
  scale_x_log10(breaks = SAMPLE_SIZES) +
  scale_colour_manual(values = mcol) + scale_fill_manual(values = mcol) +
  labs(title = "Performance scaling across the data limited range",
       subtitle = "Discrimination on interpolation tests rises with sample size and separates methods",
       x = "training sample size (log scale)", y = "test AUC",
       colour = "method", fill = "method") +
  theme_pub
savef(f1, "fig1_sample_size_scaling.png")

# ---- Figure 2  spatial structure and variance partition (RQ2) -------
vp <- as.data.table(A$variance_partition)
vp[, grp := factor(grp,
     levels = c("Residual", "cell_id", "sp_id", "dataset", "dataset:target"),
     labels = c("residual (replicate + method)", "design cell",
                "species", "dataset", "survey by target"))]
f2a <- ggplot(vp, aes(x = "", y = prop, fill = grp)) +
  geom_col(width = 0.7) +
  scale_fill_viridis_d(option = "mako", begin = 0.2, end = 0.9) +
  scale_y_continuous(labels = scales::percent) +
  coord_flip() +
  labs(title = "Where performance variation comes from",
       subtitle = "Variance partition of a hierarchical model of test AUC",
       x = NULL, y = "share of variance", fill = "source") +
  theme_pub + theme(legend.position = "bottom")

spe <- d[, .(AUC = mean(AUC)), by = .(spatial, target)]
spe[, spatial := factor(spatial, c(FALSE, TRUE),
                        c("non spatial", "spatial"))]
f2b <- ggplot(spe, aes(spatial, AUC, fill = spatial)) +
  geom_col(width = 0.6) +
  facet_wrap(~target) +
  geom_text(aes(label = sprintf("%.3f", AUC)), vjust = -0.4, size = 3.2) +
  scale_fill_viridis_d(option = "cividis", begin = 0.3, end = 0.8,
                       guide = "none") +
  coord_cartesian(ylim = c(0.5, 0.85)) +
  labs(title = "Spatial structure helps interpolation, not extrapolation",
       x = NULL, y = "mean test AUC") +
  theme_pub
f2 <- f2a / f2b + plot_layout(heights = c(1, 1.3))
savef(f2, "fig2_spatial_and_variance.png", w = 8, h = 7)

# ---- Figure 3  coverage by sample size (RQ3) ------------------------
cn <- d[grepl("marginal", target),
        .(AUC = mean(AUC)), by = .(coverage, n)]
f3 <- ggplot(cn, aes(factor(n), coverage, fill = AUC)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.2f", AUC)), size = 3.3,
            colour = "white") +
  scale_fill_viridis(option = "viridis", name = "AUC") +
  labs(title = "Spatial coverage versus sample size on extrapolation skill",
       subtitle = "Broader coverage lifts marginal AUC at every sample size",
       x = "training sample size", y = "training spatial coverage") +
  theme_pub
savef(f3, "fig3_coverage_by_sample_size.png", w = 8, h = 4.6)

# ---- Figure 4  validation optimism (RQ4, headline) ------------------
opt <- as.data.table(A$optimism)
opt[, method := factor(method, levels = opt[order(optimism), method])]
opt[, kind := factor(spatial, c(TRUE, FALSE),
                     c("spatial method", "non spatial method"))]
f4 <- ggplot(opt, aes(optimism, method, fill = kind)) +
  geom_col(width = 0.7) +
  geom_vline(xintercept = 0, colour = "grey50") +
  geom_errorbarh(aes(xmin = optimism - opt_sd / sqrt(n_pairs),
                     xmax = optimism + opt_sd / sqrt(n_pairs)),
                 height = 0.25, colour = "grey25") +
  scale_fill_manual(values = c("spatial method" = "#0072B2",
                               "non spatial method" = "#E69F00"),
                    name = NULL) +
  labs(title = "Validation optimism is largest for spatial methods",
       subtitle = "AUC lost moving from interpolation to spatial extrapolation, per method",
       x = "optimism  (conditional minus marginal AUC)", y = NULL) +
  theme_pub
savef(f4, "fig4_validation_optimism.png", w = 8, h = 4.8)

# ---- Figure 5  cross metric correlation and reliability -------------
Cm <- H$correlation_pearson
cd <- as.data.table(as.table(Cm)); setnames(cd, c("m1", "m2", "r"))
lev <- H$metrics
cd[, m1 := factor(m1, lev)][, m2 := factor(m2, rev(lev))]
f5a <- ggplot(cd, aes(m1, m2, fill = r)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f", r)), size = 3) +
  scale_fill_gradient2(low = "#b2182b", mid = "white", high = "#2166ac",
                       midpoint = 0, limits = c(-1, 1), name = "r") +
  labs(title = "Cross metric correlation", x = NULL, y = NULL) +
  theme_pub + theme(axis.text.x = element_text(angle = 45, hjust = 1))

rel <- as.data.table(H$reliability)
rel[, pair := sprintf("%s \u2192 %s", from, to)]
# a translation is classified on two criteria, the residual spread and the
# larger of its two divergences (across inference target, across survey).
# Both are shown, so the class of every translation is explicable from the
# figure rather than only from the table.
rel[, binding_div := pmax(context_divergence, survey_divergence, na.rm = TRUE)]
rel[, reliability := factor(reliability,
      levels = c("exact (analytic)", "reliable", "context dependent",
                 "unreliable"))]
rel <- rel[order(reliability, -pearson)]
rel[, pair := factor(pair, levels = rev(pair))]
relcol <- c("exact (analytic)" = "#0072B2", "reliable" = "#56B4E9",
            "context dependent" = "#E69F00", "unreliable" = "#D55E00")

L <- melt(rel[, .(pair, reliability, binding_shift, rel_resid, binding_div)],
          id.vars = c("pair", "reliability", "binding_shift"),
          measure.vars = c("rel_resid", "binding_div"),
          variable.name = "criterion", value.name = "xval")
L[, criterion := factor(criterion, levels = c("rel_resid", "binding_div"),
    labels = c("relative residual spread", "binding divergence"))]
thr <- data.table(
  criterion = factor(rep(levels(L$criterion), each = 2),
                     levels = levels(L$criterion)),
  xint = c(0.10, 0.18, 0.05, 0.10))

f5b <- ggplot(L, aes(xval, pair, colour = reliability, shape = binding_shift)) +
  geom_vline(data = thr, aes(xintercept = xint), linetype = 2,
             colour = "grey60", linewidth = 0.35) +
  geom_point(size = 2.7) +
  facet_wrap(~ criterion, scales = "free_x") +
  scale_x_log10() +
  scale_colour_manual(values = relcol, name = "translation") +
  scale_shape_manual(values = c("inference target" = 16, "survey" = 17),
                     name = "binding factor") +
  labs(x = NULL, y = NULL) +
  theme_pub
f5 <- f5a + f5b + plot_layout(widths = c(1, 1.75))
savef(f5, "fig5_metric_harmonisation.png", w = 13.5, h = 5.4)

# ---- Figure 6  calibration and overfitting (failure modes) ----------
cal <- as.data.table(A$calibration)
cal[, method := factor(method, levels = METHODS)]
f6a <- ggplot(d, aes(method, calib_slope, fill = method)) +
  geom_hline(yintercept = 1, colour = "grey40", linetype = 2) +
  geom_boxplot(outlier.size = 0.4, alpha = 0.85) +
  scale_fill_manual(values = mcol, guide = "none") +
  coord_cartesian(ylim = c(0, 3)) +
  labs(title = "Calibration slope by method",
       subtitle = "One is perfect; below one is over confident, above one under confident",
       x = NULL, y = "calibration slope") +
  theme_pub + theme(axis.text.x = element_text(angle = 30, hjust = 1))

ofm <- as.data.table(A$overfitting_method)
ofm[, method := factor(method, levels = ofm[order(gap_AUC), method])]
f6b <- ggplot(ofm, aes(gap_AUC, method, fill = method)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = mcol, guide = "none") +
  labs(title = "Optimism within the training set",
       subtitle = "Train minus test AUC gap; larger means more overfitting",
       x = "train minus test AUC", y = NULL) +
  theme_pub
f6 <- f6a + f6b + plot_layout(widths = c(1.2, 1))
savef(f6, "fig6_calibration_overfitting.png", w = 12, h = 5)

# ---- Figure 7  split failures by species rarity ---------------------
if (!is.null(A$split_failures) && nrow(A$split_failures)) {
  sf <- as.data.table(A$split_failures)
  sf[, band := factor(band, levels = c("rare", "intermediate", "common"))]
  sf2 <- sf[, .(n_fail = sum(n_fail)), by = .(band, n)]
  f7 <- ggplot(sf2, aes(factor(n), n_fail, fill = band)) +
    geom_col(position = "dodge", width = 0.75) +
    scale_fill_viridis_d(option = "inferno", begin = 0.25, end = 0.75,
                         name = "prevalence band") +
    labs(title = "Where subsampling breaks down under data limitation",
         subtitle = "Degenerate training splits concentrate in rare species at small sample sizes",
         x = "training sample size", y = "number of degenerate splits") +
    theme_pub
  savef(f7, "fig7_split_failures.png", w = 8, h = 4.6)
}

# ---------------------------------------------------------------------
# fig8  heterogeneity across surveys
# (a) validation optimism by survey and method class, against the pooled
#     values that the headline result reports
# (b) spatial advantage by survey and inference target, showing that the
#     sign of the advantage is not fixed across surveys
# Colour carries one meaning only, the method class, and matches fig4.
# Panel (b) encodes the target by facet, not by colour, so no second
# colour vocabulary is introduced.
# ---------------------------------------------------------------------
if (!is.null(A$optimism_by_dataset) && !is.null(A$spatial_effect_by_dataset)) {
  ord <- A$optimism_by_dataset_pooled[order(-optimism), dataset]
  cls_col <- c("spatial method" = "#0072B2", "non spatial method" = "#E69F00")

  o <- as.data.table(A$optimism_by_dataset)
  o[, cls := ifelse(spatial, "spatial method", "non spatial method")]
  o[, dataset := factor(dataset, levels = rev(ord))]
  # The two method classes share a survey row and their intervals overlap
  # (on NEUS the estimates differ by 0.007), so they are offset vertically
  # within the row. A deterministic offset is used in preference to random
  # jitter: an estimate carrying a confidence interval must not have its
  # plotted position perturbed at random, and the offset is reproducible.
  o[, ypos := as.numeric(dataset) + ifelse(spatial, 0.17, -0.17)]
  ref_sp   <- A$optimism[spatial == TRUE,  weighted.mean(optimism, n_pairs)]
  ref_nsp  <- A$optimism[spatial == FALSE, weighted.mean(optimism, n_pairs)]

  f8a <- ggplot(o, aes(optimism, ypos, colour = cls)) +
    geom_vline(xintercept = ref_sp,  colour = "#0072B2",
               linetype = "22", linewidth = 0.4) +
    geom_vline(xintercept = ref_nsp, colour = "#E69F00",
               linetype = "22", linewidth = 0.4) +
    geom_errorbarh(aes(xmin = optimism - 1.96 * opt_se,
                       xmax = optimism + 1.96 * opt_se),
                   height = 0.11, linewidth = 0.45) +
    geom_point(size = 2.7) +
    scale_colour_manual(values = cls_col, name = NULL) +
    scale_x_continuous(limits = c(0, NA)) +
    scale_y_continuous(breaks = seq_along(levels(o$dataset)),
                       labels = levels(o$dataset),
                       expand = expansion(add = 0.5)) +
    labs(x = "validation optimism (conditional minus marginal AUC)",
         y = NULL) +
    theme_pub

  s <- as.data.table(A$spatial_effect_by_dataset)
  s[, dataset := factor(dataset, levels = rev(ord))]
  s[, target := factor(target, levels = c("conditional", "marginal"),
                       labels = c("conditional (interpolation)",
                                  "marginal (extrapolation)"))]
  f8b <- ggplot(s, aes(spatial_advantage, dataset)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey45",
               linewidth = 0.4) +
    geom_col(fill = "#5B7183", width = 0.62) +
    facet_wrap(~ target) +
    labs(x = "spatial advantage (spatial minus non spatial AUC)",
         y = NULL) +
    theme_pub

  f8 <- (f8a / f8b) +
    plot_layout(heights = c(1, 1.05), guides = "collect") +
    plot_annotation(tag_levels = "a") &
    theme(legend.position = "bottom", legend.direction = "horizontal")
  savef(f8, "fig8_survey_heterogeneity.png", w = 8.6, h = 7.4)
}


