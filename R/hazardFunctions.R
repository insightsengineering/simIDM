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
  assert_positive_number(h01)
  assert_positive_number(h02)
  assert_positive_number(h12)

  h012 <- h12 - h01 - h02
  ((h12 - h02) * (h01 + h02) - h01 * h12 * exp(-h012 * t)) / ((h12 - h02) - h01 * exp(-h012 * t))
}

#' Helper Function for `avgHRExpOS()`
#'
#' It is an integrand of the form OS hazard function with intensities h01, h02, h12
#' at time point t multiplied with a weighted product of the two OS Survival functions
#' at t (one for intensities h0 and one for h1).
#'
#' @param x (`numeric`)\cr variable of integration.
#' @param h01 (positive `number`)\cr transition hazard for 0 to 1 transition.
#' @param h02 (positive `number`)\cr transition hazard for 0 to 2 transition.
#' @param h12 (positive `number`)\cr transition hazard for 1 to 2 transition.
#' @param h0 (`list`)\cr transition parameters for the first treatment group.
#' @param h1 (`list`)\cr transition parameters for the second treatment group.
#' @param alpha (`number`)\cr weight parameter, see [avgHRExpOS()].
#' @return This returns the value of the integrand used to calculate
#' the average hazard ratio for constant transition hazards, see [avgHRExpOS()].
#' @export
#'
#' @examples
#' h0 <- list(h01 = 0.18, h02 = 0.06, h12 = 0.17)
#' h1 <- list(h01 = 0.23, h02 = 0.07, h12 = 0.19)
#' avgHRIntegExpOS(x = 5, h01 = 0.2, h02 = 0.5, h12 = 0.7, h0 = h0, h1 = h1, alpha = 0.5)
avgHRIntegExpOS <- function(x, h01, h02, h12, h0, h1, alpha) {
  assert_positive_number(alpha)

  weightedSurv <- (ExpSurvOS(t = x, h01 = h0$h01, h02 = h0$h02, h12 = h0$h12) *
    ExpSurvOS(t = x, h01 = h1$h01, h02 = h1$h02, h12 = h1$h12))^alpha
  ExpHazOS(x, h01, h02, h12) * weightedSurv
}

#' Average OS Hazard Ratio from Constant Transition Hazards
#'
#' @param transitionByArm (`list`)\cr transition parameters for each treatment group.
#'   See [exponential_transition()], [piecewise_exponential()] and [weibull_transition()] for details.
#' @param alpha (`number`)\cr assigns weights to time points, where values higher than 0.5 assign greater weight
#' to earlier times and values lower than 0.5 assign greater weight to later times.
#' @param upper (`number`)\cr upper (time) limit of integration.
#'
#' @return This returns the value of the average hazard ratio.
#' @export
#'
#' @examples
#' transitionTrt <- exponential_transition(h01 = 0.18, h02 = 0.06, h12 = 0.17)
#' transitionCtl <- exponential_transition(h01 = 0.23, h02 = 0.07, h12 = 0.19)
#' transitionList <- list(transitionCtl, transitionTrt)
#' avgHRExpOS(transitionByArm = transitionList, alpha = 0.5, upper = 100)
avgHRExpOS <- function(transitionByArm, alpha = 0.5, upper = Inf) {
  assert_list(transitionByArm, len = 2)
  assert_class(transitionByArm[[1]], "TransitionParameters")
  assert_class(transitionByArm[[2]], "TransitionParameters")
  assert_positive_number(alpha)
  assert_positive_number(upper)

  h0 <- transitionByArm[[1]]$hazard
  h1 <- transitionByArm[[2]]$hazard

  num <- stats::integrate(avgHRIntegExpOS,
    lower = 0, upper = upper, h01 = h1$h01, h02 = h1$h02, h12 = h1$h12,
    h0 = h0, h1 = h1, alpha = alpha
  )$value
  denom <- stats::integrate(avgHRIntegExpOS,
    lower = 0, upper = upper, h01 = h0$h01, h02 = h0$h02, h12 = h0$h12,
    h0 = h0, h1 = h1, alpha = alpha
  )$value
  num / denom
}
