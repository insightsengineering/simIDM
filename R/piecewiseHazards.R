
#' Piecewise Constant Hazard Values
#'
#' This returns piecewise constant hazard values at specified time points.
#'
#' @param haz (`numeric`)\cr piecewise constant input hazard.
#' @param pw (`numeric`)\cr time intervals for the piecewise constant hazard.
#' @param x (`numeric`)\cr time-points.
#'
#' @return hazard values at input time-points.
#' @export
#'
#' @examples
#' getPCWHazard(c(1, 1.2, 1.4), c(0, 2, 3), c(1, 4, 6))
getPCWHazard <- function(haz, pw, x) {
  hazVal <- sapply(x, function(jj) {
    y <- NULL
    # find interval and corresponding hazard value for time x[jj]
    for (ii in 1:length(haz)) {
      if (jj >= pw[ii]) {
        y <- haz[ii]
      }
    }
    hazVal <- y
    return(hazVal)
  })
  return(hazVal)
}



#' Sum of two piecewise constant hazards
#'
#' Sum of two piecewise constant hazards.
#'
#' @param haz1 (`numeric`)\cr first summand (piecewise constant hazard).
#' @param haz2 (`numeric`)\cr second summand (piecewise constant hazard).
#' @param pw1 (`numeric`)\cr time intervals of first summand.
#' @param pw2  (`numeric`)\cr time intervals of second summand.
#'
#' @return constant hazards of the sum of two piecewise constant hazards.
#' @export
#'
#' @examples getSumPCW(c(1.2, 0.3, 0.6), c(1.2, 0.7, 1), c(0, 8, 9), c(0, 1, 4))
getSumPCW <- function(haz1, haz2, pw1, pw2) {
  # get all cutpoints for the intervals
  cuts_sum <- unique(sort(c(pw1, pw2)))
  haz_sum <- NULL
  ## get sum of hazards for all intervals
  for (i in 1:(length(cuts_sum))) {
    haz_sum[i] <- getPCWHazard(haz1, pw1, cuts_sum[i]) +
      getPCWHazard(haz2, pw2, cuts_sum[i])
  }
  return(haz_sum)
}
