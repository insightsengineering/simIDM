#' Event Times Distributed as Sum of Weibull
#'
#' This returns event times with a distribution resulting from the sum of two Weibull distributed random variables
#' using the inversion method.
#'
#' @param U (`numeric`)\cr uniformly distributed random variables.
#' @param haz1 (positive `number`)\cr  first summand (constant hazard).
#' @param haz2 (positive `number`)\cr  second summand (constant hazard).
#' @param p1 (positive `number`)\cr rate parameter of Weibull distribution for `haz1`.
#' @param p2 (positive `number`)\cr rate parameter of Weibull distribution for `haz2`.
#' @param entry (`numeric`)\cr the entry times in the current state.
#'
#' @return This returns a vector with event times.
#' @export
#'
#' @examples
#' getWaitTimeSum(U = c(0.4, 0.5), haz1 = 0.8, haz2 = 1, p1 = 1.1, p2 = 1.5, entry = c(0, 0))
getWaitTimeSum <- function(U, haz1, haz2, p1, p2, entry) {
  assert_numeric(U, lower = 0, upper = 1, any.missing = FALSE, all.missing = FALSE)
  assert_positive_number(haz1, zero_ok = TRUE)
  assert_positive_number(haz2, zero_ok = TRUE)
  assert_positive_number(p1)
  assert_positive_number(p2)
  assert_numeric(entry, lower = 0, len = length(U), any.missing = FALSE, all.missing = FALSE)

  N <- length(U)
  # Exponential distributed survival times.
  if (p1 == 1 && p2 == 1) {
    return(-log(1 - U) / (haz1 + haz2))
  }
  # Weibull distributed survival times.
  temp <- function(x, y, t0) {
    return(haz1 * (x + t0)^p1 - haz1 * t0^p1 - haz2 * t0^p2 + haz2 * (x + t0)^p2 + y)
  }
  stime <- NULL
  i <- 1
  while (length(stime) < N) {
    u <- U[i]
    t0 <- entry[i]
    if (temp(0, log(1 - u), t0) * temp(10000, log(1 - u), t0) < 0) {
      res <- stats::uniroot(temp, c(0, 10000), tol = 10^-16, y = log(1 - u), t0 = t0)
      stime[i] <- res$root
      i <- i + 1
      if (res$root == 0) {
        warning("uniroot: accuracy (tol=10^-16) is not enough.")
      }
    } else {
      warning("uniroot: Values at interval endpoints (interval=c(0,10000)) are not of opposite signs. \n")
    }
  }
  stime
}
