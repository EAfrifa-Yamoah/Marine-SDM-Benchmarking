\name{metric_correlations}
\alias{metric_correlations}
\title{Correlations between performance metrics in the benchmark}
\description{
Return the correlation matrix of the six performance metrics over all
benchmark fits.
}
\usage{
metric_correlations(method = c("pearson", "spearman"))
}
\arguments{
  \item{method}{\code{"pearson"} (default) or \code{"spearman"}.}
}
\value{
A 6 by 6 numeric matrix with rows and columns named \code{AUC}, \code{TSS},
\code{Spearman}, \code{R2}, \code{RMSE} and \code{Brier}.
}
\seealso{\code{\link{harmonise_metric}}}
\examples{
round(metric_correlations("spearman"), 2)
}
