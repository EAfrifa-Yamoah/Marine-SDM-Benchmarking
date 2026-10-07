#' Harmonise a species distribution model performance metric
#'
#' Translate a value on one performance metric to the expected value on
#' another, using empirical mappings learned from a controlled marine
#' species distribution model benchmark. The mapping can be conditioned
#' on the inference target, because the relationship between metrics
#' shifts between interpolation (conditional) and spatial extrapolation
#' (marginal) settings.
#'
#' A translation is certified reliable only when the mapping is stable in
#' two respects: across the inference target, and across the survey it was
#' estimated from. The second check is new in version 0.2.0. Some mappings
#' that look stable across inference targets move appreciably between
#' surveys, and treating those as portable would repeat the pooling error
#' the benchmark was built to expose. The returned \code{binding_shift}
#' column names whichever of the two is the larger source of instability.
#'
#' @param value Numeric value(s) of the source metric.
#' @param from Source metric name, one of AUC, TSS, Spearman, R2, RMSE,
#'   Brier.
#' @param to Target metric name, same set as \code{from}.
#' @param context One of "overall", "conditional" or "marginal".
#' @param prevalence Prevalence of the evaluation set on which the source
#'   value was computed. Required for an exact translation into or out of
#'   R2, by the identity R2 = 1 - Brier / (prevalence (1 - prevalence)).
#'   RMSE to Brier and Brier to RMSE are always exact, because
#'   Brier = RMSE^2, and do not use it. If omitted, R2 to RMSE returns the
#'   learned mapping, marginalised over the benchmark prevalence
#'   distribution, with a warning; RMSE to R2, Brier to R2 and R2 to Brier
#'   stop with an error asking for the prevalence.
#' @return A data frame with one row per input value and the columns value,
#'   estimate, sd (an approximate one standard deviation band), reliability,
#'   binding_shift (the larger source of instability), kind (identity,
#'   exact, exact given prevalence, exact given score distribution or
#'   empirical) and class_stability (share of
#'   species bootstrap replicates that returned the same class). Values
#'   outside the range of the source metric in the benchmark are mapped to
#'   the nearest end of that range.
#' @examples
#' harmonise_metric(0.75, from = "AUC", to = "TSS")
#' harmonise_metric(c(0.6, 0.8), "AUC", "TSS", context = "marginal")
#' @export
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

  # Brier = RMSE^2 holds for every fit, whatever the prevalence
  if (all(c(from, to) %in% c("RMSE", "Brier"))) {
    est <- if (from == "RMSE") value^2 else sqrt(pmax(value, 0))
    return(data.frame(value = value, estimate = est, sd = 0,
                      reliability = "exact", binding_shift = NA_character_,
                      kind = "exact", class_stability = NA_real_))
  }

  err <- c("R2", "RMSE", "Brier")
  if (from %in% err && to %in% err) {
    if (is.null(prevalence)) {
      if (is.null(ht$maps[[paste(from, to, sep = "__")]]))
        stop("the ", from, " to ", to, " translation needs the evaluation prevalence; ",
             "supply prevalence = <prevalence of the evaluation set>, because no learned ",
             "mapping is packaged for this direction", call. = FALSE)
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
    stop("no learned mapping for the ", from, " to ", to, " translation; ",
         "available directions: ", paste(names(ht$maps), collapse = ", "),
         call. = FALSE)

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

#' Reliability table for all learned metric translations
#' @return A data frame classifying each translation.
#' @export
harmonisation_reliability <- function() get_tables()$reliability

#' Cross metric correlation matrix from the benchmark
#' @param method "pearson" or "spearman".
#' @return A correlation matrix over the performance metrics.
#' @export
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

