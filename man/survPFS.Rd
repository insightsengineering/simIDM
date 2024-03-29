% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/corPFSOS.R
\name{survPFS}
\alias{survPFS}
\alias{survPFS.ExponentialTransition}
\alias{survPFS.WeibullTransition}
\alias{survPFS.PWCTransition}
\title{PFS Survival Function for Different Transition Models}
\usage{
survPFS(transition, t)

\method{survPFS}{ExponentialTransition}(transition, t)

\method{survPFS}{WeibullTransition}(transition, t)

\method{survPFS}{PWCTransition}(transition, t)
}
\arguments{
\item{transition}{(\code{TransitionParameters})\cr
see \code{\link[=exponential_transition]{exponential_transition()}}, \code{\link[=weibull_transition]{weibull_transition()}} or \code{\link[=piecewise_exponential]{piecewise_exponential()}} for details.}

\item{t}{(\code{numeric})\cr time at which the value of the PFS survival function is to be computed.}
}
\value{
The value of the survival function for the specified transition and time.
}
\description{
PFS Survival Function for Different Transition Models
}
\section{Methods (by class)}{
\itemize{
\item \code{survPFS(ExponentialTransition)}: Survival Function for an exponential transition model.

\item \code{survPFS(WeibullTransition)}: Survival Function for a Weibull transition model.

\item \code{survPFS(PWCTransition)}: Survival Function for a piecewise constant transition model.

}}
\examples{
transition <- exponential_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6)
survPFS(transition, 0.4)
transition <- exponential_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6)
survPFS(transition, 0.4)
transition <- weibull_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6, p01 = 2, p02 = 2.5, p12 = 3)
survPFS(transition, 0.4)
transition <- piecewise_exponential(
  h01 = c(1, 1, 1), h02 = c(1.5, 0.5, 1), h12 = c(1, 1, 1),
  pw01 = c(0, 3, 8), pw02 = c(0, 6, 7), pw12 = c(0, 8, 9)
)
survPFS(transition, 0.4)
}
