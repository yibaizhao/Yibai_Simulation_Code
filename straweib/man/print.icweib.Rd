\name{print.icweib}
\alias{print.icweib}
\title{
Print icweib object
}
\description{
This function prints the summary of the fitting results from icweib. 
}
\usage{
\method{print}{icweib}(x, digits = 3, ...)
}

\arguments{
  \item{x}{
output returned by icweib function. 
}
  \item{digits}{
  digits to print.
}
  \item{\dots}{
  other arguments to be passed from print function.
}
}


\author{
Xiangdong Gu and Raji Balasubramanian
}

\examples{
data(tooth24)
fit <- icweib(L = left, R = right, data = tooth24, strata = dmf, covariates = ~sex)
fit
}

