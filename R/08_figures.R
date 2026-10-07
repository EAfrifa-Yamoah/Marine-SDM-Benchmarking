# =====================================================================
# 08_figures.R
# Result figures of the benchmark, written as PNG to figures/.
# Colour follows R/palette.R: hue identifies a method and nothing else;
# every other distinction is carried by position, facet, shape, line
# type or a direct label in neutral greys. No figure carries a title;
# captions belong in the manuscript text.
# =====================================================================

suppressMessages({
  library(data.table); library(ggplot2); library(viridis); library(patchwork)
})
source(file.path("R", "config.R"))
source(file.path("R", "palette.R"))

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
mcol <- METHOD_COLOURS[METHODS]
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
  labs(x = "training sample size (log scale)", y = "test AUC",
       colour = "method", fill = "method") +
  theme_pub
savef(f1, "fig1_sample_size_scaling.png")

# ---- Figure 2  spatial structure and variance partition (RQ2) -------
# The variance sources are named on the axis and the shares printed, so
# no fill colour is needed to tell them apart.
vp <- as.data.table(A$variance_partition)
vp[, grp := factor(grp,
     levels = rev(c("Residual", "sp_id", "dataset", "cell_id", "dataset:target")),
     labels = rev(c("residual (replicate and method)", "species", "survey",
                    "design cell", "survey by target")))]
f2a <- ggplot(vp, aes(prop, grp)) +
  geom_col(width = 0.62, fill = INK_MID) +
  geom_text(aes(label = sprintf("%.1f%%", 100 * prop)), hjust = -0.15,
            size = 3.3) +
  scale_x_continuous(labels = scales::percent,
                     expand = expansion(mult = c(0, 0.12))) +
  labs(x = "share of variance in test AUC", y = NULL) +
  theme_pub

spe <- d[, .(AUC = mean(AUC)), by = .(spatial, target)]
spe[, spatial := factor(spatial, c(FALSE, TRUE),
                        c("non-spatial methods", "spatial methods"))]
f2b <- ggplot(spe, aes(spatial, AUC)) +
  geom_col(width = 0.6, fill = INK_MID) +
  facet_wrap(~target) +
  geom_text(aes(label = sprintf("%.3f", AUC)), vjust = -0.4, size = 3.2) +
  coord_cartesian(ylim = c(0.5, 0.85)) +
  labs(x = NULL, y = "mean test AUC") +
  theme_pub
f2 <- f2a / f2b + plot_layout(heights = c(1, 1.3)) +
  plot_annotation(tag_levels = "a")
savef(f2, "fig2_spatial_and_variance.png", w = 8, h = 7)

# ---- Figure 3  coverage by sample size (RQ3) ------------------------
cn <- d[grepl("marginal", target),
        .(AUC = mean(AUC)), by = .(coverage, n)]
f3 <- ggplot(cn, aes(factor(n), coverage, fill = AUC)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.3f", AUC)), size = 3.3,
            colour = "white") +
  scale_fill_viridis(option = "viridis", name = "mean\nAUC") +
  labs(x = "training sample size", y = "training spatial coverage") +
  theme_pub
savef(f3, "fig3_coverage_by_sample_size.png", w = 8, h = 4.6)

# ---- Figure 4  validation optimism (RQ4, headline) ------------------
# Bars take the method colours of Figure 1; the method class is shown by
# the panel each method sits in, not by a second colour vocabulary.
# intervals from species means: fits within a species are not independent,
# so the standard error of a method mean is taken over its 60 species means
pw <- dcast(d, dataset + species + n + coverage + rep + method + spatial ~ target, value.var = "AUC")
setnames(pw, c("conditional (interpolation)", "marginal (extrapolation)"), c("cond", "marg"))
pw <- pw[!is.na(cond) & !is.na(marg)][, o := cond - marg]
sp_means <- pw[, .(o = mean(o)), by = .(method, spatial, dataset, species)]
opt <- sp_means[, .(optimism = mean(o), se = sd(o) / sqrt(.N)), by = .(method, spatial)]
opt <- merge(opt[, .(method, spatial, se)], as.data.table(A$optimism)[, .(method = as.character(method), optimism)], by = "method")
opt[, method := factor(method, levels = opt[order(optimism), method])]
opt[, kind := factor(spatial, c(TRUE, FALSE),
                     c("spatial methods", "non-spatial methods"))]
f4 <- ggplot(opt, aes(optimism, method, fill = method)) +
  geom_col(width = 0.7) +
  geom_vline(xintercept = 0, colour = "grey50") +
  geom_errorbarh(aes(xmin = optimism - 1.96 * se, xmax = optimism + 1.96 * se),
                 height = 0.25, colour = "grey25") +
  facet_grid(kind ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(values = mcol, guide = "none") +
  labs(x = "optimism (conditional minus marginal AUC)", y = NULL) +
  theme_pub + theme(strip.text.y = element_text(angle = 0, hjust = 0))
savef(f4, "fig4_validation_optimism.png", w = 8, h = 4.8)

# ---- Figure 5  cross metric correlation and reliability -------------
Cm <- H$correlation_pearson
cd <- as.data.table(as.table(Cm)); setnames(cd, c("m1", "m2", "r"))
lev <- H$metrics
cd[, m1 := factor(m1, lev)][, m2 := factor(m2, rev(lev))]
f5a <- ggplot(cd, aes(m1, m2, fill = r)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f", r),
                colour = r > 0.55), size = 3, show.legend = FALSE) +
  scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = "white")) +
  scale_fill_viridis(option = "viridis", limits = c(-0.3, 1), name = "r") +
  labs(x = NULL, y = NULL) +
  theme_pub + theme(axis.text.x = element_text(angle = 45, hjust = 1))

rel <- as.data.table(H$reliability)
# a translation is classified on two criteria, the residual spread and the
# larger of its two divergences (across inference target, across survey).
# Both are shown, so the class of every translation is explicable from the
# figure rather than only from the table. The class is carried by shape;
# an asterisk marks translations whose binding divergence is the survey.
rel[, binding_div := pmax(context_divergence, survey_divergence, na.rm = TRUE)]
rel[, pair := sprintf("%s → %s%s", from, to,
                      ifelse(binding_shift == "survey", " *", ""))]
rel[, reliability := factor(reliability,
      levels = c("exact (analytic)", "reliable", "context dependent",
                 "unreliable"))]
rel <- rel[order(reliability, -pearson)]
rel[, pair := factor(pair, levels = rev(pair))]

L <- melt(rel[, .(pair, reliability, rel_resid, binding_div)],
          id.vars = c("pair", "reliability"),
          measure.vars = c("rel_resid", "binding_div"),
          variable.name = "criterion", value.name = "xval")
L[, criterion := factor(criterion, levels = c("rel_resid", "binding_div"),
    labels = c("relative residual spread", "binding divergence"))]
thr <- data.table(
  criterion = factor(rep(levels(L$criterion), each = 2),
                     levels = levels(L$criterion)),
  xint = c(0.10, 0.18, 0.05, 0.10))

f5b <- ggplot(L, aes(xval, pair, shape = reliability)) +
  geom_vline(data = thr, aes(xintercept = xint), linetype = 2,
             colour = "grey60", linewidth = 0.35) +
  geom_point(size = 2.7, colour = INK, stroke = 0.9) +
  facet_wrap(~ criterion, scales = "free_x") +
  scale_x_log10() +
  scale_shape_manual(values = RELIABILITY_SHAPES, name = "translation") +
  labs(x = NULL, y = NULL) +
  theme_pub
f5 <- f5a + f5b + plot_layout(widths = c(1, 1.75))
savef(f5, "fig5_metric_harmonisation.png", w = 13.5, h = 5.4)

# ---- Figure 6  calibration and overfitting (failure modes) ----------
f6a <- ggplot(d, aes(method, calib_slope, fill = method)) +
  geom_hline(yintercept = 1, colour = "grey40", linetype = 2) +
  geom_boxplot(outlier.shape = NA, alpha = 0.85) +
  scale_fill_manual(values = mcol, guide = "none") +
  coord_cartesian(ylim = c(0, 3)) +
  labs(x = NULL, y = "calibration slope") +
  theme_pub + theme(axis.text.x = element_text(angle = 30, hjust = 1))

ofm <- as.data.table(A$overfitting_method)
ofm[, method := factor(method, levels = ofm[order(gap_AUC), method])]
f6b <- ggplot(ofm, aes(gap_AUC, method, fill = method)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = mcol, guide = "none") +
  labs(x = "train minus test AUC", y = NULL) +
  theme_pub
f6 <- f6a + f6b + plot_layout(widths = c(1.2, 1)) +
  plot_annotation(tag_levels = "a")
savef(f6, "fig6_calibration_overfitting.png", w = 12, h = 5)

# ---- Figure 7  split failures by species rarity ---------------------
# One panel per prevalence band, each on its own scale, so the band is
# named by the panel label and no fill colour is needed.
if (!is.null(A$split_failures) && nrow(A$split_failures)) {
  sf <- as.data.table(A$split_failures)
  sf2 <- sf[, .(n_fail = sum(n_fail)), by = .(band, n)]
  sf2[, band := factor(band, levels = c("rare", "intermediate", "common"),
                       labels = c("rare species", "intermediate species",
                                  "common species"))]
  f7 <- ggplot(sf2, aes(factor(n), n_fail)) +
    geom_col(width = 0.7, fill = INK_MID) +
    geom_text(aes(label = format(n_fail, big.mark = ",")), vjust = -0.4,
              size = 3.1) +
    facet_wrap(~ band, scales = "free_y") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(x = "training sample size", y = "number of degenerate splits") +
    theme_pub
  savef(f7, "fig7_split_failures.png", w = 8, h = 4.2)
}

# ---------------------------------------------------------------------
# fig8  heterogeneity across surveys
# (a) validation optimism by survey and method class, against the pooled
#     values that the headline result reports; class by shape, pooled
#     references by line type
# (b) spatial advantage by survey and inference target
# ---------------------------------------------------------------------
if (!is.null(A$optimism_by_dataset) && !is.null(A$spatial_effect_by_dataset)) {
  ord <- A$optimism_by_dataset_pooled[order(-optimism), dataset]

  o <- as.data.table(A$optimism_by_dataset)
  cl <- pw[, .(o = mean(o)), by = .(dataset, spatial, species)][, .(se_cl = sd(o) / sqrt(.N)), by = .(dataset, spatial)]
  o <- merge(o, cl, by = c("dataset", "spatial"))
  o[, cls := factor(ifelse(spatial, "spatial method", "non-spatial method"),
                    levels = names(CLASS_SHAPES))]
  o[, dataset := factor(dataset, levels = rev(ord))]
  # The two method classes share a survey row and their intervals can
  # overlap, so they are offset vertically within the row by a fixed,
  # reproducible amount rather than by random jitter.
  o[, ypos := as.numeric(dataset) + ifelse(spatial, 0.17, -0.17)]
  ref <- data.table(
    cls = factor(names(CLASS_SHAPES), levels = names(CLASS_SHAPES)),
    x = c(A$optimism[spatial == TRUE,  weighted.mean(optimism, n_pairs)],
          A$optimism[spatial == FALSE, weighted.mean(optimism, n_pairs)]))

  f8a <- ggplot(o, aes(optimism, ypos, shape = cls)) +
    geom_vline(data = ref, aes(xintercept = x, linetype = cls),
               colour = INK_MID, linewidth = 0.45) +
    geom_errorbarh(aes(xmin = optimism - 1.96 * se_cl,
                       xmax = optimism + 1.96 * se_cl),
                   height = 0.11, linewidth = 0.45, colour = INK) +
    geom_point(size = 2.7, colour = INK, fill = "white", stroke = 0.9) +
    scale_shape_manual(values = CLASS_SHAPES, name = NULL) +
    scale_linetype_manual(values = c("spatial method" = "22",
                                     "non-spatial method" = "12"),
                          name = NULL) +
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
    geom_col(fill = INK_MID, width = 0.62) +
    facet_wrap(~ target) +
    labs(x = "spatial advantage (spatial minus non-spatial AUC)",
         y = NULL) +
    theme_pub

  f8 <- (f8a / f8b) +
    plot_layout(heights = c(1, 1.05), guides = "collect") +
    plot_annotation(tag_levels = "a") &
    theme(legend.position = "bottom", legend.direction = "horizontal",
          legend.key.width = unit(1.6, "lines"))
  savef(f8, "fig8_survey_heterogeneity.png", w = 8.6, h = 7.4)
}
