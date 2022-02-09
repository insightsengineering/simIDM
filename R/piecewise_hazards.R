## for piecewise-constant hazards: function to get hazard value at time x

#' Piecewise Constant Hazard Values
#'
#' This calculates piecewise constant hazard values ...
#'
#' @param haz (`numeric`)\cr input hazards.
#' @param pw bla.
#' @param x bla.
#'
#' @return Bli.
#' @export
#'
#' @examples
#' getPCWHazard(1:5, 2:3, 2)
getPCWHazard <- function(haz, pw, x) {
  hazVal <- sapply(x, function(jj) {
    y <- NULL
    # find interval and correponding hazard value for time x[jj]
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
