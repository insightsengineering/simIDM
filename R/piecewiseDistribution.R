#' Piecewise Exponentially Distributed Event Times
#'
#' This returns event times with a distribution resulting from piece-wise constant hazards
#' using the inversion method.
#'
#' @param U (`numeric`)\cr uniformly distributed random variables.
#' @param haz (`numeric`)\cr piecewise constant hazard.
#' @param pw (`numeric`)\cr time intervals for the piecewise constant hazard.
#' @param t_0 (`numeric`)\cr the starting times.
#'
#' @return This returns a vector with event times.
#' @export
#'
#' @examples
#' getPCWDistr(U = runif(3), haz = c(1.1, 0.5, 0.4), pw = c(0, 7, 10), t_0 = c(0, 1, 4.2))
getPCWDistr <- function(U, haz, pw, t_0) {
  assert_numeric(U, lower = 0, upper = 1, any.missing = FALSE, all.missing = FALSE)
  assert_numeric(haz, lower = 0, any.missing = FALSE, all.missing = FALSE)
  assert_intervals(pw, length(haz))
  assert_numeric(t_0, lower = 0, len = length(U), any.missing = FALSE, all.missing = FALSE)

  N <- length(U)
  t1 <- rep(NA, N) # Initialize event times.
  n2 <- length(haz)
  for (kk in seq_len(N)) {
    # Shift if t_0 !=0.
    cuts_temp <- c(t_0[kk], pw[pw > t_0[kk]])
    cuts_temp <- cuts_temp - t_0[kk]
    # Number of cut-points after shift.
    n <- length(cuts_temp)
    # Number of hazards to remove for shift (t_0 > times).
    remov <- n2 - n
    haz_temp <- haz[(remov + 1):n2]
    LogU <- log(1 - U[kk])

    t1[kk] <- if (n != 1) {
      PCWInversionMethod(haz = haz_temp, pw = cuts_temp, LogU = LogU)
    } else if (n == 1) { # I.e. exponential distribution.
      -LogU / haz_temp
    }
  }
  t1
}

#' Single Piecewise Exponentially Distributed Event Time
#'
#' This returns an event time with a distribution resulting from piece-wise constant hazards
#' using the inversion method.
#'
#' @param pw (`numeric`)\cr time intervals for the piecewise constant hazard.
#' @param haz (`numeric`)\cr piecewise constant hazard.
#' @param LogU (`numeric`)\cr transformed uniformly distributed random variables (log(1-U)).
#'
#' @return This returns one single event time.
#' @export
#'
#' @examples
#' PCWInversionMethod(haz = c(1.1, 0.5, 0.4), pw = c(0, 7, 10), LogU = log(1 - runif(1)))
PCWInversionMethod <- function(haz, pw, LogU) {
  assert_numeric(haz, lower = 0, any.missing = FALSE, all.missing = FALSE)
  assert_intervals(pw, length(haz))
  assert_number(LogU)

  n <- length(pw)
  t1 <- 0
  # Determine sum of alpha*time-interval for all i.
  dt <- pw[2:n] - pw[1:(n - 1)]

  # Helping matrix.
  tempMatrix <- matrix(0, nrow = n, ncol = n - 1)
  tempMatrix[lower.tri(tempMatrix)] <- 1

  sumA <- -as.vector(tempMatrix %*% (haz[1:(n - 1)] * dt))
  # Find the appropriate time interval.
  rev_index <- findInterval(LogU, rev(sumA), rightmost.closed = FALSE, left.open = TRUE, all.inside = FALSE)
  # Use reverse index.
  index <- n - rev_index
  # Calculate event time.
  pw[index] + (sumA[index] - LogU) / haz[index]
}
