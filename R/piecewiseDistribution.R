
#' Piecewise exponentially distributed event times
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
#' getPCWDistr(runif(3), c(1.1, 0.5, 0.4), c(0, 7, 10), c(0, 1, 4.2))
getPCWDistr <- function(U, haz, pw, t_0) {
  N <- length(U)
  t1 <- rep(NA, N) # Initialize event times.
  n2 <- length(haz)
  for (kk in seq_len(N)) {
    # Shift if t_0 !=0.
    cuts_temp <- c(t_0[kk], pw[pw > t_0[kk]])
    cuts_temp <- cuts_temp - t_0[kk]
    # Number of cutpoints after shift.
    n <- length(cuts_temp)
    # Number of hazards to remove for shift (t_0 > times).
    remov <- n2 - n
    haz_temp <- haz[(remov + 1):n2]
    LogU <- log(1 - U[kk])

    if (n != 1) {
      ## Determine sum of alpha*time-interval for all i.
      dt <- cuts_temp[2:n] - cuts_temp[1:(n - 1)]

      # Helping matrix.
      tempMatrix <- matrix(0, nrow = n, ncol = n - 1)
      tempMatrix[lower.tri(tempMatrix)] <- 1

      sumA <- -as.vector(tempMatrix %*% (haz_temp[1:(n - 1)] * dt))

      # Find the appropriate time interval.
      for (i in 1:n) {
        if (i != n) {
          t1[kk] <- ifelse(sumA[i] >= LogU & LogU > sumA[i + 1],
            cuts_temp[i] + (sumA[i] - LogU) / haz_temp[i],
            t1[kk]
          )
        } else {
          t1[kk] <- ifelse(LogU <= sumA[i],
            cuts_temp[i] + (sumA[i] - LogU) / haz_temp[i],
            t1[kk]
          )
        }
      }
    } else if (n == 1) { # I.e. exponential distribution.
      t1[kk] <- -LogU / haz_temp
    }
  }
  return(t1)
}
