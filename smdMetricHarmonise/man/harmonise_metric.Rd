\name{harmonise_metric}
\alias{harmonise_metric}
\title{Harmonise a species distribution model performance metric}
\description{
Translate a value on one performance metric to the expected value on
another, using mappings learned from a controlled benchmark of species
distribution models fitted to four marine bottom trawl surveys.
}
\usage{
harmonise_metric(value, from, to, context = "overall", prevalence = NULL)
}
\arguments{
  \item{value}{Numeric value(s) of the source metric.}
  \item{from}{Source metric: one of \code{"AUC"}, \code{"TSS"},
    \code{"Spearman"}, \code{"R2"}, \code{"RMSE"} or \code{"Brier"}.}
  \item{to}{Target metric, from the same set.}
  \item{context}{Inference target under which the source value was
    obtained: \code{"conditional"} (interpolation, for example random
    cross-validation), \code{"marginal"} (spatial extrapolation, for example
    a held-out spatial block) or \code{"overall"} (both pooled).}
  \item{prevalence}{Prevalence of the evaluation set on which the source
    value was computed, a number or vector strictly between 0 and 1. With
    it, translations into or out of \code{"R2"} are exact, by the identity
    R2 = 1 - Brier / (prevalence (1 - prevalence)). Without it, R2 to RMSE
    returns the learned mapping with a warning, and RMSE to R2, Brier to R2
    and R2 to Brier stop with an error. RMSE to Brier and Brier to RMSE are
    always exact, because Brier = RMSE^2, and ignore this argument.}
}
\details{
RMSE to Brier and Brier to RMSE are returned exactly by Brier = RMSE^2.
Translations into or out of R2 among the three error and explained
variation metrics are returned exactly when \code{prevalence} is supplied.
Learned mappings are used for AUC to TSS, TSS to AUC, AUC to Spearman,
Spearman to AUC, TSS to Spearman, AUC to R2, AUC to RMSE, AUC to Brier and,
without \code{prevalence}, R2 to RMSE. The packaged tables also hold learned
mappings for RMSE to Brier and Brier to RMSE, which are kept for reference
and not used. Any other pair of different metrics stops with an error
listing the available directions.

Each learned translation carries a reliability class. A translation is
\code{"reliable"} when its residual standard deviation is below 0.10 of the
range of the target metric and the larger of its two shifts (between
inference targets and between surveys) is below 0.05, \code{"context
dependent"} when these are below 0.18 and 0.10, and \code{"unreliable"}
otherwise. In \code{\link{harmonisation_reliability}} the learned Brier to
RMSE and RMSE to Brier mappings are marked \code{"exact (analytic)"}; an
exact translation returned by this function is marked \code{"exact"}. An
unreliable translation returns its estimate with a warning. \code{binding_shift} names the larger of the two
shifts, and \code{class_stability} is the share of 100 cluster bootstrap
replicates over species within survey that returned the same class.

A value outside the range of the source metric observed in the benchmark is
mapped to the translation at the nearest end of that range; see the package
vignette for a way to flag such values.
}
\value{
A data frame with one row per input value and the columns \code{value},
\code{estimate}, \code{sd} (an approximate one standard deviation band, the
residual standard deviation of the learned mapping, or 0 for exact
translations), \code{reliability}, \code{binding_shift}, \code{kind}
(\code{"identity"}, \code{"exact"}, \code{"exact given prevalence"},
\code{"exact given score distribution"} or \code{"empirical"}) and
\code{class_stability}.
}
\seealso{
\code{\link{harmonisation_reliability}}, \code{\link{metric_correlations}},
\code{\link{harmonise_tables}}, and the tutorial in
\code{vignette("smdMetricHarmonise")}.
}
\examples{
harmonise_metric(0.75, from = "AUC", to = "TSS")
harmonise_metric(c(0.6, 0.8), "AUC", "TSS", context = "marginal")

# RMSE and Brier are related exactly
harmonise_metric(0.30, from = "RMSE", to = "Brier")

# exact translation when the evaluation prevalence is known
harmonise_metric(0.20, from = "R2", to = "Brier", prevalence = 0.25)

# an unreliable translation returns its estimate with a warning
suppressWarnings(harmonise_metric(0.75, from = "AUC", to = "R2"))
}
