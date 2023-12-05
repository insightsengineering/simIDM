#' Piecewise Constant Hazard Values
#'
#' This returns piecewise constant hazard values at specified time points.
#'
#' @param haz (`numeric`)\cr piecewise constant input hazard.
#' @param pw (`numeric`)\cr time intervals for the piecewise constant hazard.
#' @param x (`numeric`)\cr time-points.
#'
#' @return Hazard values at input time-points.
#' @export
#'
#' @examples
#' getPWCHazard(c(1, 1.2, 1.4), c(0, 2, 3), c(1, 4, 6))
getPWCHazard <- function(haz, pw, x) { # nolint
  sapply(x, function(jj) {
    y <- NULL
    # Find interval and corresponding hazard value for time x[jj].
    for (ii in seq_along(haz)) {
      if (jj >= pw[ii]) {
        y <- haz[ii]
      }
    }
    y
  })
}

#' Sum of Two Piecewise Constant Hazards
#'
#' This returns the sum of two piecewise constant hazards per interval.
#'
#' @param haz1 (`numeric`)\cr first summand (piecewise constant hazard).
#' @param haz2 (`numeric`)\cr second summand (piecewise constant hazard).
#' @param pw1 (`numeric`)\cr time intervals of first summand.
#' @param pw2  (`numeric`)\cr time intervals of second summand.
#'
#' @return List with elements `hazards` and `intervals` for the sum of two piecewise constant hazards.
#' @export
#'
#' @examples
#' getSumPCW(c(1.2, 0.3, 0.6), c(1.2, 0.7, 1), c(0, 8, 9), c(0, 1, 4))
getSumPCW <- function(haz1, haz2, pw1, pw2) {
  # Get all cutpoints for the intervals.
  cuts_sum <- unique(sort(c(pw1, pw2)))
  haz_sum <- NULL
  # Get sum of hazards for all intervals.
  for (i in seq_along(cuts_sum)) {
    haz_sum[i] <- getPWCHazard(haz1, pw1, cuts_sum[i]) +
      getPWCHazard(haz2, pw2, cuts_sum[i])
  }
  list(
    hazards = haz_sum,
    intervals = cuts_sum
  )
}
