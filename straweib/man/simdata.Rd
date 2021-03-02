\name{simdata}
\alias{simdata}
\docType{data}
\title{
Simulated data with mixed types of events
}
\description{
This simulated data contains event times that are left censored, right censored, interval censored, or non-censored (observed event). The data is generated from a stratified Weibull distribution model in which each stratum is assumed to have an independent stratum-specific shape parameter. In addition, the regression coefficients corresponding to the vector of explanatory variables excluding the stratum indicator are assumed to be constant across strata.
}
\usage{data(simdata)}
\format{
  A data frame with 298 observations on the following 6 variables:
  \describe{
    \item{\code{ID}}{subject id}
    \item{\code{strata}}{strata}
    \item{\code{cov1}}{a continuous covariate}
    \item{\code{cov2}}{a continuous covariate}
    \item{\code{left}}{left endpoint of censoring interval}
    \item{\code{right}}{right endpoint of censoring interval}
  }
}
\references{
see \code{\link{icweib}} for details on how to set \emph{L} and \emph{R} for different types of events.
}
\examples{
data(simdata)
}

