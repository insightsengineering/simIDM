% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/corPFSOS.R
\name{survPFSOS}
\alias{survPFSOS}
\title{Survival Function of the Product PFS*OS for Different Transition Models}
\usage{
survPFSOS(t, transition)
}
\arguments{
\item{t}{(\code{numeric})\cr time at which the value of the PFS*OS survival function is to be computed.}

\item{transition}{(\code{TransitionParameters})\cr
see \code{\link[=exponential_transition]{exponential_transition()}}, \code{\link[=weibull_transition]{weibull_transition()}} or \code{\link[=piecewise_exponential]{piecewise_exponential()}} for details.}
}
\value{
This returns the value of PFS*OS survival function at time t.
}
\description{
Survival Function of the Product PFS*OS for Different Transition Models
}
\examples{
transition <- exponential_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6)
survPFSOS(0.4, transition)
}
