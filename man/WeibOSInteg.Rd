% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/survivalFunctions.R
\name{WeibOSInteg}
\alias{WeibOSInteg}
\title{Helper Function for \code{WeibSurvOS()}}
\usage{
WeibOSInteg(x, t, h01, h02, h12, p01, p02, p12)
}
\arguments{
\item{x}{(\code{numeric})\cr  variable of integration.}

\item{t}{(\code{numeric})\cr  study time-points.}

\item{h01}{(positive \code{number})\cr transition hazard for 0 to 1 transition.}

\item{h02}{(positive \code{number})\cr transition hazard for 0 to 2 transition.}

\item{h12}{(positive \code{number})\cr transition hazard for 1 to 2 transition.}

\item{p01}{(positive \code{number})\cr rate parameter of Weibull distribution for \code{h01}.}

\item{p02}{(positive \code{number})\cr rate parameter of Weibull distribution for \code{h02}.}

\item{p12}{(positive \code{number})\cr rate parameter of Weibull distribution for \code{h12}.}
}
\value{
Numeric results of the integrand used to calculate
the OS survival function for Weibull transition hazards, see  \code{WeibSurvOS()}.
}
\description{
Helper Function for \code{WeibSurvOS()}
}
\keyword{internal}
