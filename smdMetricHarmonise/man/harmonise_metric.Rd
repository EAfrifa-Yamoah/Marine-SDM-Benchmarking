\name{harmonise_metric}
\alias{harmonise_metric}
\title{Harmonise a species distribution model performance metric}
\description{Translate a value on one performance metric to the
expected value on another, using empirical mappings learned from a
marine species distribution model benchmark.}
\usage{harmonise_metric(value, from, to, context = "overall")}
\arguments{
  \item{value}{Numeric value(s) of the source metric.}
  \item{from}{Source metric: AUC, TSS, Spearman, R2, RMSE or Brier.}
  \item{to}{Target metric, same set as from.}
  \item{context}{One of overall, conditional or marginal.}}
\value{A data frame with the translated estimate, an approximate
standard deviation band and a reliability class.}
\examples{harmonise_metric(0.75, from = "AUC", to = "TSS")}
