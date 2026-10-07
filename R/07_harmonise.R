# =====================================================================
# 07_harmonise.R
# Empirical cross metric harmonisation protocol (research question 5).
#
# Studies report species distribution model skill in different currencies:
# some in discrimination (AUC, TSS), some in explained variation (R2),
# some in error (RMSE, Brier). A synthesis needs to translate between
# them. This script learns the empirical mapping between every metric
# pair from the benchmark fits, measures how tight and how context stable
# each mapping is, and classifies each translation as reliable, context
# dependent or unreliable. The learned tables and a lookup function are
# assembled into the smdMetricHarmonise package.
# =====================================================================

suppressMessages({ library(data.table); library(mgcv) })
source(file.path("R", "config.R"))

RES <- PATHS$results
PKG <- file.path(dirname(RES), "smdMetricHarmonise")
d <- as.data.table(readRDS(file.path(RES, "fit_metrics.rds")))

METS <- c("AUC", "TSS", "Spearman", "R2", "RMSE", "Brier")
d[, target := factor(target, levels = c("conditional", "marginal"))]

# ---- correlation structure ------------------------------------------
Cpear <- cor(d[, ..METS], use = "pairwise.complete.obs")
Cspear <- cor(d[, ..METS], use = "pairwise.complete.obs", method = "spearman")
fwrite(as.data.table(round(Cpear, 3), keep.rownames = "metric"),
       file.path(RES, "table_metric_correlations.csv"))

# ---- learn pairwise mappings ----------------------------------------
# For an ordered pair (from -> to) a monotone friendly smooth of the
# target metric on the source metric is fitted, overall and separately by
# inference target. The prediction grid, the residual spread and the
# divergence between the conditional and marginal mappings are stored.
grid_n <- 50L
learn_pair <- function(from, to) {
  sub <- d[is.finite(d[[from]]) & is.finite(d[[to]])]
  if (nrow(sub) < 50) return(NULL)
  rng <- quantile(sub[[from]], c(0.02, 0.98), names = FALSE)
  xg <- seq(rng[1], rng[2], length.out = grid_n)

  fit_on <- function(dat) {
    f <- as.formula(sprintf("%s ~ s(%s, k=5)", to, from))
    m <- tryCatch(gam(f, data = dat, method = "REML"),
                  error = function(e) NULL)
    if (is.null(m)) return(NULL)
    nd <- setNames(data.frame(xg), from)
    list(pred = as.numeric(predict(m, nd)),
         resid_sd = sd(residuals(m)))
  }
  overall <- fit_on(sub)
  cond <- fit_on(sub[target == "conditional"])
  marg <- fit_on(sub[target == "marginal"])
  if (is.null(overall)) return(NULL)

  # context divergence: mean absolute gap between conditional and marginal
  # predicted target values over the shared source range
  ctx_div <- if (!is.null(cond) && !is.null(marg))
    mean(abs(cond$pred - marg$pred)) else NA_real_

  # survey divergence: the same idea applied to the survey of origin. Each
  # survey mapping is predicted only over its own support, and the
  # divergence is the mean pairwise absolute difference among the survey
  # curves across the range where every survey has data. For two groups
  # this definition reduces exactly to the context divergence above, so the
  # two quantities are on one scale and share the same thresholds. A
  # translation that is stable between inference targets can still move
  # between surveys, and pooling the surveys hides that.
  dss <- unique(sub$dataset)
  P <- vapply(dss, function(ds) {
    sd_ <- sub[dataset == ds]
    if (nrow(sd_) < 50) return(rep(NA_real_, grid_n))
    q  <- quantile(sd_[[from]], c(0.02, 0.98), names = FALSE)
    fo <- fit_on(sd_)
    if (is.null(fo)) return(rep(NA_real_, grid_n))
    p <- fo$pred
    p[xg < q[1] | xg > q[2]] <- NA_real_      # no extrapolation off support
    p
  }, numeric(grid_n))
  ok <- stats::complete.cases(P)
  sur_div <- if (sum(ok) >= 5)
    mean(apply(P[ok, , drop = FALSE], 1,
               function(r) mean(stats::dist(r)))) else NA_real_

  list(from = from, to = to, x = xg,
       y_overall = overall$pred,
       y_conditional = if (!is.null(cond)) cond$pred else NA,
       y_marginal = if (!is.null(marg)) marg$pred else NA,
       resid_sd = overall$resid_sd,
       context_divergence = ctx_div,
       survey_divergence = sur_div,
       survey_support = sum(ok),
       pearson = unname(Cpear[from, to]),
       spearman = unname(Cspear[from, to]),
       to_range = diff(range(sub[[to]], na.rm = TRUE)))
}

pairs <- list(
  c("AUC", "TSS"), c("TSS", "AUC"),
  c("AUC", "Spearman"), c("Spearman", "AUC"),
  c("TSS", "Spearman"),
  c("AUC", "R2"), c("AUC", "RMSE"), c("AUC", "Brier"),
  c("RMSE", "Brier"), c("Brier", "RMSE"),
  c("R2", "RMSE"))
maps <- Filter(Negate(is.null), lapply(pairs, function(p) learn_pair(p[1], p[2])))
names(maps) <- vapply(maps, function(m) paste(m$from, m$to, sep = "__"),
                      character(1))

# ---- classify reliability -------------------------------------------
# A translation is reliable when the residual spread is a small fraction of
# the target metric range and the mapping barely moves between validation
# settings. Two settings are checked, the inference target and the survey of
# origin, and the binding one is the larger of the two divergences. A
# mapping that is stable across inference targets but moves between surveys
# is not portable, and reporting it as reliable would repeat in the protocol
# the same pooling error the benchmark exposes in Section 3.7.
classify <- function(m) {
  rel_sd <- m$resid_sd / m$to_range
  ctx <- m$context_divergence
  sur <- m$survey_divergence
  if (is.na(ctx)) ctx <- Inf
  if (is.na(sur)) sur <- Inf          # unassessable stability is not stability
  shift <- max(ctx, sur)
  if (rel_sd < 0.10 && shift < 0.05) "reliable"
  else if (rel_sd < 0.18 && shift < 0.10) "context dependent"
  else "unreliable"
}
reliab <- rbindlist(lapply(maps, function(m) data.table(
  from = m$from, to = m$to, pearson = round(m$pearson, 3),
  spearman = round(m$spearman, 3),
  resid_sd = round(m$resid_sd, 4),
  rel_resid = round(m$resid_sd / m$to_range, 3),
  context_divergence = round(m$context_divergence, 4),
  survey_divergence = round(m$survey_divergence, 4),
  binding_shift = ifelse(m$survey_divergence > m$context_divergence,
                         "survey", "inference target"),
  reliability = classify(m))))

# analytic note: Brier equals RMSE squared for probability predictions, so
# that translation is exact rather than empirical
reliab[from == "RMSE" & to == "Brier", reliability := "exact (analytic)"]
reliab[from == "Brier" & to == "RMSE", reliability := "exact (analytic)"]
setorder(reliab, reliability, -pearson)
fwrite(reliab, file.path(RES, "table_harmonisation_reliability.csv"))

log_file <- file.path(PATHS$logs, "harmonise.log"); cat("", file = log_file)
say <- function(...) { m <- sprintf(...); cat(m, "\n"); cat(m, "\n", file = log_file, append = TRUE) }
say("=== harmonisation over %d fits ===", nrow(d))
say("reliable translations: %s",
    paste(reliab[reliability == "reliable",
                 sprintf("%s->%s", from, to)], collapse = ", "))
say("unreliable translations: %s",
    paste(reliab[reliability == "unreliable",
                 sprintf("%s->%s", from, to)], collapse = ", "))
say("translations whose binding instability is the survey, not the inference target: %s",
    paste(reliab[binding_shift == "survey",
                 sprintf("%s->%s (survey %.3f vs context %.3f)",
                         from, to, survey_divergence, context_divergence)],
          collapse = "; "))

harmonise_tables <- list(
  maps = maps, reliability = reliab,
  correlation_pearson = Cpear, correlation_spearman = Cspear,
  metrics = METS, built = as.character(Sys.Date()),
  n_fits = nrow(d), n_surveys = length(unique(d$dataset)),
  theory = if (file.exists(file.path(RES, "table_harmonisation_theory.csv")))
    as.data.frame(fread(file.path(RES, "table_harmonisation_theory.csv"))) else NULL,
  version = "0.3.0")
saveRDS(harmonise_tables, file.path(RES, "harmonise_tables.rds"))

# =====================================================================
# assemble the smdMetricHarmonise package
# =====================================================================
dir.create(file.path(PKG, "R"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(PKG, "data"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(PKG, "man"), showWarnings = FALSE, recursive = TRUE)

# package data
save(harmonise_tables, file = file.path(PKG, "data", "harmonise_tables.rda"),
     compress = "xz")

# The package's R code, help pages, DESCRIPTION, NAMESPACE and tutorial vignette are
# maintained as source files in smdMetricHarmonise/; this script refreshes only its data.
say("smdMetricHarmonise data written to %s", file.path(PKG, "data", "harmonise_tables.rda"))
