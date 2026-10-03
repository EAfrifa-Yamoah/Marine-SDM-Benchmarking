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

# DESCRIPTION
writeLines(c(
  "Package: smdMetricHarmonise",
  "Type: Package",
  "Title: Harmonising Species Distribution Model Performance Metrics",
  "Version: 0.3.0",
  "Authors@R: person('Eben', 'Afrifa-Yamoah', role = c('aut', 'cre'),",
  "    email = 'e.afrifayamoah@ecu.edu.au')",
  paste("Description: Empirical translation between species distribution",
        "model performance metrics (AUC, TSS, Spearman, R2, RMSE, Brier)",
        "learned from a controlled benchmark across four marine survey",
        "datasets under data limited conditions. Provides a lookup that",
        "maps a value on one metric to the expected value on another,",
        "with a reliability classification indicating when a translation",
        "is trustworthy, context dependent or should be avoided. A",
        "translation is certified reliable only when the mapping is stable",
        "both across the inference target, that is interpolation against",
        "spatial extrapolation, and across the survey of origin. The",
        "survey stability check was added in 0.2.0. Version 0.3.0 returns",
        "translations among the coefficient of determination, root mean",
        "squared error and Brier score by identity when the evaluation",
        "prevalence is supplied, since R2 = 1 - Brier / (pi (1 - pi)) holds",
        "exactly, and warns when it is not; it also reports whether a",
        "learned class is stable under a cluster bootstrap over species."),
  "License: MIT + file LICENSE",
  "Encoding: UTF-8",
  "LazyData: true",
  "Depends: R (>= 4.0)",
  "RoxygenNote: 7.2.3"),
  file.path(PKG, "DESCRIPTION"))

writeLines(c(
  "YEAR: 2026",
  "COPYRIGHT HOLDER: Eben Afrifa-Yamoah"),
  file.path(PKG, "LICENSE"))

# NAMESPACE
writeLines(c(
  "export(harmonise_metric)",
  "export(harmonisation_reliability)",
  "export(metric_correlations)"),
  file.path(PKG, "NAMESPACE"))

# R source
pkg_r <- '#\' Harmonise a species distribution model performance metric
#\'
#\' Translate a value on one performance metric to the expected value on
#\' another, using empirical mappings learned from a controlled marine
#\' species distribution model benchmark. The mapping can be conditioned
#\' on the inference target, because the relationship between metrics
#\' shifts between interpolation (conditional) and spatial extrapolation
#\' (marginal) settings.
#\'
#\' A translation is certified reliable only when the mapping is stable in
#\' two respects: across the inference target, and across the survey it was
#\' estimated from. The second check is new in version 0.2.0. Some mappings
#\' that look stable across inference targets move appreciably between
#\' surveys, and treating those as portable would repeat the pooling error
#\' the benchmark was built to expose. The returned \\code{binding_shift}
#\' column names whichever of the two is the larger source of instability.
#\'
#\' @param value Numeric value(s) of the source metric.
#\' @param from Source metric name, one of AUC, TSS, Spearman, R2, RMSE,
#\'   Brier.
#\' @param to Target metric name, same set as \\code{from}.
#\' @param context One of "overall", "conditional" or "marginal".
#\' @param prevalence Prevalence of the evaluation set on which the source
#\'   value was computed. Required for an exact translation among R2, RMSE
#\'   and Brier, which are related by the identity R2 = 1 - Brier /
#\'   (prevalence (1 - prevalence)). If omitted the learned mapping,
#\'   marginalised over the benchmark prevalence distribution, is returned
#\'   with a warning.
#\' @return A data frame with the input value, the translated estimate, an
#\'   approximate one standard deviation band, the reliability class and the
#\'   binding source of instability.
#\' @examples
#\' harmonise_metric(0.75, from = "AUC", to = "TSS")
#\' harmonise_metric(c(0.6, 0.8), "AUC", "TSS", context = "marginal")
#\' @export
harmonise_metric <- function(value, from, to, context = "overall", prevalence = NULL) {
  stopifnot(is.numeric(value))
  ht <- get_tables()
  from <- match.arg(from, ht$metrics)
  to   <- match.arg(to,   ht$metrics)
  context <- match.arg(context, c("overall", "conditional", "marginal"))
  if (identical(from, to))
    return(data.frame(value = value, estimate = value, sd = 0,
                      reliability = "identity", binding_shift = NA_character_,
                      kind = "identity", class_stability = NA_real_))

  err <- c("R2", "RMSE", "Brier")
  if (from %in% err && to %in% err) {
    if (is.null(prevalence)) {
      warning("the ", from, " to ", to, " translation is exact given the evaluation ",
              "prevalence, which was not supplied; the learned mapping marginalised ",
              "over the benchmark prevalence distribution is returned instead",
              call. = FALSE)
    } else {
      stopifnot(is.numeric(prevalence), all(prevalence > 0 & prevalence < 1))
      v  <- prevalence * (1 - prevalence)
      bs <- switch(from, Brier = value, RMSE = value^2, R2 = (1 - value) * v)
      est <- switch(to, Brier = bs, RMSE = sqrt(pmax(bs, 0)), R2 = 1 - bs / v)
      return(data.frame(value = value, estimate = est, sd = 0,
                        reliability = "exact", binding_shift = NA_character_,
                        kind = "exact given prevalence", class_stability = NA_real_))
    }
  }

  key <- paste(from, to, sep = "__")
  m <- ht$maps[[key]]
  if (is.null(m))
    stop("no learned mapping for ", from, " -> ", to,
         "; available: ", paste(names(ht$maps), collapse = ", "))

  yv <- switch(context,
               overall = m$y_overall,
               conditional = m$y_conditional,
               marginal = m$y_marginal)
  if (length(yv) == 1L && is.na(yv[1])) yv <- m$y_overall
  est <- stats::approx(m$x, yv, xout = value, rule = 2)$y

  row <- ht$reliability[ht$reliability$from == from &
                          ht$reliability$to == to, ]
  rc <- if (nrow(row)) as.character(row$reliability[1]) else NA_character_
  bs <- if (nrow(row)) as.character(row$binding_shift[1]) else NA_character_
  th <- ht$theory
  trow <- if (!is.null(th)) th[th$from == from & th$to == to, ] else NULL
  kd <- if (!is.null(trow) && nrow(trow)) as.character(trow$kind[1]) else "empirical"
  st <- if (!is.null(trow) && nrow(trow) && "stab_range" %in% names(trow)) trow$stab_range[1] else NA_real_

  if (identical(rc, "unreliable"))
    warning("the ", from, " to ", to, " translation is classified unreliable",
            if (!is.na(bs)) paste0(" (destabilised chiefly by the ", bs, ")"),
            "; the estimate is returned but should not be used for synthesis",
            call. = FALSE)

  data.frame(value = value, estimate = est,
             sd = m$resid_sd,
             reliability = rc, binding_shift = bs,
             kind = kd, class_stability = st)
}

#\' Reliability table for all learned metric translations
#\' @return A data frame classifying each translation.
#\' @export
harmonisation_reliability <- function() get_tables()$reliability

#\' Cross metric correlation matrix from the benchmark
#\' @param method "pearson" or "spearman".
#\' @return A correlation matrix over the performance metrics.
#\' @export
metric_correlations <- function(method = c("pearson", "spearman")) {
  method <- match.arg(method)
  ht <- get_tables()
  if (method == "pearson") ht$correlation_pearson else ht$correlation_spearman
}

# internal: load the packaged tables, falling back for interactive use
get_tables <- function() {
  e <- new.env()
  utils::data("harmonise_tables", package = "smdMetricHarmonise",
              envir = e)
  e$harmonise_tables
}
'
writeLines(pkg_r, file.path(PKG, "R", "harmonise.R"))

# minimal manual page for the main function
writeLines(c(
  "\\name{harmonise_metric}",
  "\\alias{harmonise_metric}",
  "\\title{Harmonise a species distribution model performance metric}",
  "\\description{Translate a value on one performance metric to the",
  "expected value on another, using empirical mappings learned from a",
  "marine species distribution model benchmark.}",
  "\\usage{harmonise_metric(value, from, to, context = \"overall\")}",
  "\\arguments{",
  "  \\item{value}{Numeric value(s) of the source metric.}",
  "  \\item{from}{Source metric: AUC, TSS, Spearman, R2, RMSE or Brier.}",
  "  \\item{to}{Target metric, same set as from.}",
  "  \\item{context}{One of overall, conditional or marginal.}}",
  "\\value{A data frame with the translated estimate, an approximate",
  "standard deviation band and a reliability class.}",
  "\\examples{harmonise_metric(0.75, from = \"AUC\", to = \"TSS\")}"),
  file.path(PKG, "man", "harmonise_metric.Rd"))

say("smdMetricHarmonise package written to %s", PKG)
