% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hazardFunctions.R
\name{avgHRExpOS}
\alias{avgHRExpOS}
\title{Average OS Hazard Ratio from Constant Transition Hazards}
\usage{
avgHRExpOS(transitionByArm, alpha = 0.5, upper = Inf)
}
\arguments{
\item{transitionByArm}{(\code{list})\cr transition parameters for each treatment group.
See \code{\link[=exponential_transition]{exponential_transition()}}, \code{\link[=piecewise_exponential]{piecewise_exponential()}} and \code{\link[=weibull_transition]{weibull_transition()}} for details.}

\item{alpha}{(\code{number})\cr assigns weights to time points, where values higher than 0.5 assign greater weight
to earlier times and values lower than 0.5 assign greater weight to later times.}

\item{upper}{(\code{number})\cr upper (time) limit of integration.}
}
\value{
This returns the value of the average hazard ratio.
}
\description{
Average OS Hazard Ratio from Constant Transition Hazards
}
\examples{
transitionTrt <- exponential_transition(h01 = 0.18, h02 = 0.06, h12 = 0.17)
transitionCtl <- exponential_transition(h01 = 0.23, h02 = 0.07, h12 = 0.19)
transitionList <- list(transitionCtl, transitionTrt)
avgHRExpOS(transitionByArm = transitionList, alpha = 0.5, upper = 100)
}
