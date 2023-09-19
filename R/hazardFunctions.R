#' OS Hazard Function from Constant Transition Hazards
#'
#' @param t (`numeric`)\cr study time-points.
#' @param h01 (positive `number`)\cr transition hazard for 0 to 1 transition.
#' @param h02 (positive `number`)\cr transition hazard for 0 to 2 transition.
#' @param h12 (positive `number`)\cr transition hazard for 1 to 2 transition.
#'
#' @return This returns the value of the OS hazard function at time t.
#' @export
#'
#' @examples
#' ExpHazOS(c(1:5), 0.2, 1.1, 0.8)
ExpHazOS <- function(t, h01, h02, h12) {
  assert_numeric(t, lower = 0, any.missing = FALSE)
  assert_positive_number(h01, zero_ok = TRUE)
  assert_positive_number(h02, zero_ok = TRUE)
  assert_positive_number(h12, zero_ok = TRUE)

  h012 <- h12 - h01 - h02
  ((h12 - h02) * (h01 + h02) - h01 * h12 * exp(-h012 * t)) / ((h12 - h02) - h01 * exp(-h012 * t))
}
