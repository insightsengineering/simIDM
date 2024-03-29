% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/empSignificant.R
\name{logRankTest}
\alias{logRankTest}
\title{Log-Rank Test for a Single Trial}
\usage{
logRankTest(data, typeEvent = c("PFS", "OS"), critical)
}
\arguments{
\item{data}{(\code{data.frame})\cr data frame containing entry and exit times of an
illness-death model. See \code{\link[=getSimulatedData]{getSimulatedData()}} for details.}

\item{typeEvent}{(\code{string})\cr endpoint to be evaluated, possible values are \code{PFS} and \code{OS}.}

\item{critical}{(positive \code{number})\cr critical value of the log-rank test.}
}
\value{
Logical value indicating log-rank test significance.
}
\description{
This function evaluates the significance of either PFS or OS endpoints in a trial,
based on a pre-specified critical value.
}
\examples{
transition1 <- exponential_transition(h01 = 0.06, h02 = 0.3, h12 = 0.3)
transition2 <- exponential_transition(h01 = 0.1, h02 = 0.4, h12 = 0.3)
simTrial <- getClinicalTrials(
  nRep = 1, nPat = c(800, 800), seed = 1234, datType = "1rowPatient",
  transitionByArm = list(transition1, transition2), dropout = list(rate = 0.5, time = 12),
  accrual = list(param = "intensity", value = 7)
)[[1]]
logRankTest(data = simTrial, typeEvent = "OS", critical = 3.4)
}
