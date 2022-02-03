#' Exponential Transition Hazards
#'
#' This creates a list with class `TransitionParameters` containing
#' hazards, intervals and Weibull rates for exponential survival times
#' in an illness-death model.
#'
#' @param h01 (positive `number`)\cr transition hazard for 0 to 1 transition.
#' @param h02 (positive `number`)\cr transition hazard for 0 to 2 transition.
#' @param h12 (positive `number`)\cr transition hazard for 1 to 2 transition.
#'
#' @return List with elements `hazards`, `intervals`, `weibull_rates` and `family`
#'   (exponential).
#' @export
#'
#' @examples
#' exponential_transition(1, 1.6, 0.3)
exponential_transition <- function(h01, h02, h12) {
  assert_positive_number(h01, zero_ok = TRUE)
  assert_positive_number(h02, zero_ok = TRUE)
  assert_positive_number(h12, zero_ok = TRUE)
  structure(
    list(
      hazards = list(h01 = h01, h02 = h02, h12 = h12),
      intervals = list(pw01 = 0, pw02 = 0, pw12 = 0),
      weibull_rates = list(p01 = 1, p02 = 1, p12 = 1),
      family = "exponential"
    ),
    class = "TransitionParameters"
  )
}


#' todo: similar as above
#'
#' Class of transition parameters describing
#' transition hazards for piecewise exponentially distributed survival times
#' in an Illness-Death model
#'
#' @param h01 vector of constant transition hazards for 0 to 1 transition
#' @param h02 vector of constant transition hazards for 0 to 2 transition
#' @param h12 vector of constant transition hazards for 1 to 2 transition
#' @param pw01 vector of time intervals for the piecewise constant hazards `h01`
#' @param pw02 vector of time intervals for the piecewise constant hazards `h02`
#' @param pw12 vector of time intervals for the piecewise constant hazards `h12`
#'
#' @return
#' list of transition parameters for piecewise exponentially distributed
#' survival times in an Illness-Death model
#' @export
#'
#' @examples
#' piecewise_exponential(
#'   h01 = c(1, 1, 1), h02 = c(1.5, 0.5, 1), h12 = c(1, 1, 1),
#'   pw01 = c(0, 3, 8), pw02 = c(0, 6, 7), pw12 = c(0, 8, 9)
#' )
piecewise_exponential <- function(h01, h02, h12, pw01, pw02, pw12) {
  assert_numeric(h01, lower = 0, any.missing = FALSE, all.missing = FALSE)
  assert_numeric(h02, lower = 0, any.missing = FALSE, all.missing = FALSE)
  assert_numeric(h12, lower = 0, any.missing = FALSE, all.missing = FALSE)

  assert_intervals(pw01, length(h01))
  assert_intervals(pw02, length(h02))
  assert_intervals(pw12, length(h12))

  structure(
    list(
      hazards = list(h01 = h01, h02 = h02, h12 = h12),
      intervals = list(pw01 = pw01, pw02 = pw02, pw12 = pw12),
      weibull_rates = list(p01 = 1, p02 = 1, p12 = 1),
      family = "piecewise exponential"
    ),
    class = "TransitionParameters"
  )
}



#' todo: similar as above
#' Class of transition parameters describing Weibull distributed survival times
#' in an illness-death model
#'
#' @param h01 transition hazard for 0 to 1 transition
#' @param h02 transition hazard for 0 to 2 transition
#' @param h12 transition hazard for 1 to 2 transition
#' @param p01 rate parameter of Weibull distribution for `h01`
#' @param p02 rate parameter of Weibull distribution for `h02`
#' @param p12 rate parameter of Weibull distribution for `h12`
#'
#' @return
#' list of transition parameters for Weibull distributed survival times in
#'  an illness-death model
#' @export
#'
#' @examples
#' weibull_transition(h01 = 1, h02 = 1.3, h12 = 0.5, p01 = 1.2, p02 = 1.3, p12 = 0.5)
weibull_transition <- function(h01, h02, h12, p01, p02, p12) {
  assert_positive_number(h01, zero_ok = TRUE)
  assert_positive_number(h02, zero_ok = TRUE)
  assert_positive_number(h12, zero_ok = TRUE)
  assert_positive_number(p01)
  assert_positive_number(p02)
  assert_positive_number(p12)
  structure(
    list(
      hazards = list(h01 = h01, h02 = h02, h12 = h12),
      weibull_rates = list(p01 = p01, p02 = p02, p12 = p12),
      intervals = list(pw01 = 0, pw02 = 0, pw12 = 0),
      family = "Weibull"
    ),
    class = "TransitionParameters"
  )
}
