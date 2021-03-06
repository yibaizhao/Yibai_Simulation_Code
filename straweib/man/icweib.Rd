\name{icweib}
\alias{icweib}

\title{
Fit stratified Weibull regression model
}
\description{
This function fits a stratified Weibull regression model using maximum likelihood estimation. The function can incorporate right, left, interval censored outcomes in addition to fully observed (i.e. uncensored) time to event data. (see details). 

}
\usage{
icweib(L, R, data, strata = "ALL", covariates = NULL)
}
\arguments{
  \item{L}{
left endpoint of censoring interval. To indicate left censoring, set L=0.
}
  \item{R}{
right endpoint of censoring interval. To indicate right censoring, set R=Inf.
}
  \item{data}{
dataset
}
  \item{strata}{
variable for stratification. Set it to a character or numeric constant for un-stratified model.
}
  \item{covariates}{
a formula to specify explanatory variables in the proportional hazards model. The input is a right hand formula object, such as ~x1 + x2. The default is NULL, corresponding to the no covariate case.
}
}
\details{
As in the stratified Cox proportional hazards model (Collett (2003)), this model allows a baseline hazard function that is stratum-specific. However, the model assumes that the regression coefficients for all other explanatory variables (excluding the stratum indicator) are constant across strata.  Assuming a Weibull distribution for the random variable corresponding to the time to event in conjunction with the Cox proportional hazards model, the survival function can be expressed as S(t | Z) = exp(-lambda*exp(beta*Z)*t^(gamma)), where Z denotes the vector of covariates, gamma denotes the shape parameter and lambda the scale parameter.  To allow stratum-specific baseline hazard functions, we generalize the model given above by  expressing the survival function as S(t | Z, Stratum=i) = exp(-lambda_i*exp(beta*Z)*t^(gamma_i)), where i denotes the stratum, Z denotes the vector of covariates, gamma_i and lambda_i denote the shape and scale parameters for stratum i, respectively.  In particular, the model assumes that the coefficients for explanatory covariates Z (denoted by beta) are the same for all strata i. 

In the likelihood optimization, u_i=log(lambda_i) and v_i=log(gamma_i) are used as parameters to remove the parameters' range constriction. The likelihood function is optimized using optim() function. The maximum likelihood estimates are used to estimate baseline hazard ratios between two subjects (see \code{\link{HRatio}}), and survival function (see \code{\link{plot.icweib}}).

This function can accommondate different types of censored time-to-event outcomes: left censoring, right censoring, interval censoring, and non-censoring (event), by appropriately setting L and R,

\tabular{lll}{
	 L     \tab R     \tab INTERPRETATION \cr
	 \cr
	 a   \tab b   \tab interval censoring, [a, b] \cr
	 0   \tab b   \tab left censoring, [0, b] \cr
	 a   \tab Inf \tab right censoring, [a, Inf] \cr
	 a   \tab a   \tab no censoring, event time = a \cr
  } 

}
\value{
This function returns an object with class \emph{icweib}. The items in the object are,
\item{loglik}{log-likelihood functions of the full, reduced, and null models. Reduced model refers to the model that all shape parameters are same. Null model refers to the model that there is no covariate in the model.}
\item{coef}{results for estimated coefficients for explanatory variables.}
\item{weib}{estimated Weibull shape and scale parameters for each stratum.}
\item{stratatest}{results of likelihood ratio test and Wald test corresponding to the null hypothesis that all the strata specific shape parameters are equal.}
\item{cov}{covariance matrix of the parameters}
\item{ns}{information of different counts}
\item{delete}{observation numbers in the data that are deleted due to inappropriate input.}
\item{maxT}{maximum observed time in the data}
\item{q}{returned object from the \emph{optim} function for the full model.}


}

\author{
Xiangdong Gu and Raji Balasubramanian
}

\seealso{
\code{\link{HRatio}}, \code{\link{plot.icweib}}
}

\references{
Collett, D. (2003). \emph{Modelling Survival Data in Medical Research, Second Edition}, Texts in statistical science. Taylor & Francis.
}

\examples{
## Analyze tooth data
data(tooth24)   ## load data
## Stratified on dmf, and sex as explanatory variable 
fit <- icweib(L = left, R = right, data = tooth24, strata = dmf, covariates = ~sex)

## Analyze hypernephroma data
data(hyper)

## Derive left and right endpoints from time and status
hyper$left <- hyper$time
hyper$right <- ifelse(hyper$status==1, hyper$time, Inf)

## Stratified on nephrectomy, and age group as explanatory variable
fit1 <- icweib(L = left, R = right, data = hyper, strata = nephrectomy, covariates = ~factor(age)) 

}

