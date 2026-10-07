\name{harmonisation_reliability}
\alias{harmonisation_reliability}
\title{Reliability of the learned metric translations}
\description{
Return the table that classifies each of the eleven learned translations.
}
\usage{
harmonisation_reliability()
}
\value{
A data frame with one row per translation and the columns \code{from},
\code{to}, \code{pearson} and \code{spearman} (correlations between the two
metrics over the benchmark fits), \code{resid_sd} (residual standard
deviation of the learned mapping, in units of the target metric),
\code{rel_resid} (that residual as a fraction of the range of the target
metric), \code{context_divergence} (shift of the mapping between the
conditional and marginal inference targets), \code{survey_divergence}
(mean pairwise shift of the mapping between surveys), \code{binding_shift}
(the larger of the two shifts) and \code{reliability} (the class).
}
\seealso{\code{\link{harmonise_metric}}, \code{\link{harmonise_tables}}}
\examples{
rel <- harmonisation_reliability()
rel[, c("from", "to", "rel_resid", "binding_shift", "reliability")]
}
