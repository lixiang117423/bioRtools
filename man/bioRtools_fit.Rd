\name{bioRtools_fit}
\alias{bioRtools_fit}
\title{Shared result contract for bioRtools analyses}
\value{
A list-like object with class \code{bioRtools_fit}.
}
\description{
Core multivariate workflows return a common list-like object with four
standard components: \code{model}, \code{tables}, \code{plots}, and
\code{provenance}. Historical list elements remain available for backward
compatibility.
}
\details{
\itemize{
  \item \code{model}: the native model object or named list of native fits.
  \item \code{tables}: named analysis-ready data frames.
  \item \code{plots}: named ggplot objects.
  \item \code{provenance}: analysis type, engine package and version, package
  version, R version, validated input dimensions, call, and optional seed.
}
}
\seealso{
\code{\link{as_tibble.bioRtools_fit}}
}
