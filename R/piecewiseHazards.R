
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
    for (ii in seq_len(haz)) {
      if (jj >= pw[ii]) {
        y <- haz[ii]
      }
    }
    hazVal <- y
    return(hazVal)
  })
  return(hazVal)
}
