% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/survivalFunctions.R
\name{singleExpQuantOS}
\alias{singleExpQuantOS}
\title{Helper Function for Single Quantile for OS Survival Function}
\usage{
singleExpQuantOS(q, h01, h02, h12)
}
\arguments{
\item{q}{(\code{number})\cr single quantile at which to compute event time.}

\item{h01}{(\verb{numeric vector})\cr constant transition hazards for 0 to 1 transition.}

\item{h02}{(\verb{numeric vector})\cr constant transition hazards for 0 to 2 transition.}

\item{h12}{(\verb{numeric vector})\cr constant transition hazards for 1 to 2 transition.}
}
\value{
Single time t such that the OS survival function at t equals q.
}
\description{
Helper Function for Single Quantile for OS Survival Function
}
\keyword{internal}
