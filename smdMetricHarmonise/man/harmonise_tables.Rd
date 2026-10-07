\name{harmonise_tables}
\alias{harmonise_tables}
\docType{data}
\title{Translation tables learned from the benchmark}
\description{
The mappings, reliability classes and correlations used by
\code{\link{harmonise_metric}}, learned from a factorial benchmark of seven
species distribution modelling methods fitted to four FISHGLOB bottom trawl
surveys.
}
\usage{harmonise_tables}
\format{
A list with the components
\describe{
  \item{\code{maps}}{a list of eleven learned mappings named
    \code{"FROM__TO"}; each holds a grid \code{x} of 50 values of the source
    metric spanning its range in the benchmark, the fitted target values
    \code{y_overall}, \code{y_conditional} and \code{y_marginal} on that grid,
    the residual standard deviation \code{resid_sd}, the two shifts, the
    correlations and the range of the target metric \code{to_range}.}
  \item{\code{reliability}}{the table returned by
    \code{\link{harmonisation_reliability}}.}
  \item{\code{correlation_pearson}, \code{correlation_spearman}}{the
    matrices returned by \code{\link{metric_correlations}}.}
  \item{\code{metrics}}{the six metric names.}
  \item{\code{built}}{the date the tables were built.}
  \item{\code{n_fits}}{the number of benchmark fits the tables were learned
    from.}
  \item{\code{n_surveys}}{the number of surveys.}
  \item{\code{theory}}{a data frame describing all sixteen translations
    studied in the benchmark, including those returned exactly: their kind,
    residual spread under range and interdecile normalisation, target and
    survey shifts, the lower bound on the between survey residual, the class
    under both normalisation rules and its stability under the species
    bootstrap.}
  \item{\code{version}}{the package version the tables were built for.}
}
}
\source{
Afrifa-Yamoah, E. (2026). Marine SDM Benchmarking (Version 1.0.0)
[Computer software]. Zenodo. \doi{10.5281/zenodo.23202487}
}
\examples{
harmonise_tables$n_fits
names(harmonise_tables$maps)
range(harmonise_tables$maps$AUC__TSS$x)
}
\keyword{datasets}
