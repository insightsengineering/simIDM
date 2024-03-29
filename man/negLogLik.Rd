% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimateParams.R
\name{negLogLik}
\alias{negLogLik}
\title{Compute the Negative Log-Likelihood for a Given Data Set and Transition Model}
\usage{
negLogLik(transition, data)
}
\arguments{
\item{transition}{(\code{ExponentialTransition} or \code{WeibullTransition})\cr
see \code{\link[=exponential_transition]{exponential_transition()}} or \code{\link[=weibull_transition]{weibull_transition()}} for details.}

\item{data}{(\code{data.frame})\cr in the format created by \code{\link[=prepareData]{prepareData()}}.}
}
\value{
The value of the negative log-likelihood.
}
\description{
Compute the Negative Log-Likelihood for a Given Data Set and Transition Model
}
\details{
Calculates the negative log-likelihood for a given data set and transition model. It uses the hazard
and survival functions specific to the transition model.
}
\examples{
transition <- exponential_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6)
simData <- getOneClinicalTrial(
  nPat = c(30), transitionByArm = list(transition),
  dropout = list(rate = 0.8, time = 12),
  accrual = list(param = "time", value = 1)
)
negLogLik(transition, prepareData(simData))
}
