\name{plot.icweib}
\alias{plot.icweib}
\title{
Plot estimated survival function
}
\description{
This function plots the estimated survival function along with associated pointwise 95\% confidence intervals corresponding to the input strata and values of explanatory variables.  
}
\usage{
\method{plot}{icweib}(x, strata = NULL, Z = NULL, tRange = NULL, tEst = NULL, ...)
}

\arguments{
  \item{x}{
output returned by icweib function.
}
  \item{strata}{
the vector of strata for which the survival function is estimated and plotted. The default is NULL, corresponding to all strata. 
}
  \item{Z}{
the vector of values of explanatory variables for which the survival function is estimated and plotted. The order and length should match the estimated coefficients as shown in x$coef. The default is NULL, corresponding to all 0 or baseline. 
}
  \item{tRange}{
the range of time to plot. It should be in the format of c(t1, t2), which means that the range of time is t1 to t2. The default is NULL, corresponding to 0 to maximum observed time in the data. 
}
  \item{tEst}{
the vector of times at which the survival function along with associated pointwise 95\% confidence internval is estimated and output. The default is NULL, which means no estimated survival function is output.    	
}
  \item{\dots}{
arguments of \code{\link{plot}} function except \emph{col} and \emph{lty}. For example, the axis labels and title of the plot can be specified.
}
}
\details{
The survival function and associated pointwise 95\% confidence intervals are estimated for input values of time and covariates. 
}
\value{
If tEst is specified, then a dataframe of estimated survival function along with 95\% confidence interval is returned.
}
\author{
Xiangdong Gu and Raji Balasubramanian
}

\seealso{
\code{\link{icweib}} \code{\link{plot}}
}
\examples{
data(tooth24)
fit <- icweib(L = left, R = right, data = tooth24, strata = dmf, covariates = ~sex)
surv <- plot(fit, Z = 1, tRange = c(1, 7), tEst=1:7, xlab = "Time", ylab = "Survival Function", 
     main = "Estimated survival function for sex = 1 (girls)")
}

