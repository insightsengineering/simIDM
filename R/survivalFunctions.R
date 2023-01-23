#' PFS Survival Function from Constant Transition Hazards
#'
#' @inheritParams WeibSurvOS
#' @return This returns the value of PFS survival function at time t.
#' @export
#'
#' @examples
#' ExpSurvPFS(c(1:5), 0.2, 0.4)
ExpSurvPFS <- function(t, h01, h02) {
  assert_numeric(t, lower = 0)
  assert_positive_number(h01, zero_ok = TRUE)
  assert_positive_number(h02, zero_ok = TRUE)

  exp(-(h01 + h02) * t)
}

#' OS Survival Function from Constant Transition Hazards
#'
#' @inheritParams WeibSurvOS
#'
#' @return This returns the value of OS survival function at time t.
#' @export
#'
#' @examples
#' ExpSurvOS(c(1:5), 0.2, 0.4, 0.1)
ExpSurvOS <- function(t, h01, h02, h12) {
  assert_numeric(t, lower = 0)
  assert_positive_number(h01, zero_ok = TRUE)
  assert_positive_number(h02, zero_ok = TRUE)
  assert_positive_number(h12, zero_ok = TRUE)

  h012 <- h12 - h01 - h02
  ExpSurvPFS(t, h01, h02) / h012 * (h12 - h02 - h01 * exp(-h012 * t))
}

#' PFS Survival Function from Weibull Transition Hazards
#'
#' @inheritParams WeibSurvOS
#' @return This returns the value of PFS survival function at time t.
#' @export
#' @export
#'
#' @examples
#' WeibSurvPFS(c(1:5), 0.2, 0.5, 1.2, 0.9)
WeibSurvPFS <- function(t, h01, h02, p01, p02) {
  assert_numeric(t, lower = 0)
  assert_positive_number(h01, zero_ok = TRUE)
  assert_positive_number(h02, zero_ok = TRUE)
  assert_positive_number(p01)
  assert_positive_number(p02)

  exp(-h01 * t^p01 - h02 * t^p02)
}

#' Helper for Efficient Integration
#'
#' @param integrand (`function`)\cr  to be integrated.
#' @param upper (`numeric`)\cr upper limits of integration.
#' @param ...  additional arguments to be passed to `integrand`.
#'
#' @return This function returns for each upper limit the estimates of the integral.
#' @export
#'
#' @examples
#' integrand <- function(x) x^2
#' upper <- c(0, 1, 0.4, 2, 5, 2, 0.3, 0.4, 1)
#' integrateVector(integrand, upper = upper)
integrateVector <- function(integrand, upper, ...) {
  assert_true(all(upper >= 0))
  boundaries <- sort(unique(upper))
  nIntervals <- length(boundaries)
  intervals <- mapply(
    FUN = function(f, lower, upper, ...) stats::integrate(f, ..., lower, upper)$value,
    MoreArgs = list(f = integrand, ...),
    lower = c(0, boundaries[-nIntervals]),
    upper = boundaries
  )
  cumIntervals <- cumsum(intervals)
  cumIntervals[match(upper, boundaries)]
}

#' Helper Function for `WeibSurvOS()`
#'
#' @param x (`numeric`)\cr  variable of integration.
#' @inheritParams WeibSurvOS
#'
#' @return Numeric results of the integrand used to calculate
#' the OS survival function for Weibull transition hazards, see  `WeibSurvOS()`.
#' @export
#'
#' @examples
#' WeibOSInteg(1:5, 0.2, 0.5, 2.1, 1.2, 0.9, 1)
WeibOSInteg <- function(x, h01, h02, h12, p01, p02, p12) {
  x^(p01 - 1) * exp(-h01 * x^p01 - h02 * x^p02 + h12 * x^p12)
}

#' OS Survival Function from Weibull Transition Hazards
#'
#' @param t  (`numeric`)\cr  study time-points.
#' @param h01 (positive `number`)\cr transition hazard for 0 to 1 transition.
#' @param h02 (positive `number`)\cr transition hazard for 0 to 2 transition.
#' @param h12 (positive `number`)\cr transition hazard for 1 to 2 transition.
#' @param p01 (positive `number`)\cr rate parameter of Weibull distribution for `h01`.
#' @param p02 (positive `number`)\cr rate parameter of Weibull distribution for `h02`.
#' @param p12 (positive `number`)\cr rate parameter of Weibull distribution for `h12`.
#'
#' @return This returns the value of OS survival function at time t.
#' @export
#'
#' @examples WeibSurvOS(c(1:5), 0.2, 0.5, 2.1, 1.2, 0.9, 1)
WeibSurvOS <- function(t, h01, h02, h12, p01, p02, p12) {
  assert_numeric(t, lower = 0)
  assert_positive_number(h01, zero_ok = TRUE)
  assert_positive_number(h02, zero_ok = TRUE)
  assert_positive_number(h12, zero_ok = TRUE)
  assert_positive_number(p01)
  assert_positive_number(p02)
  assert_positive_number(p12)

  WeibSurvPFS(t, h01, h02, p01, p02) +
    h01 * p01 * exp(-h12 * t^p12) *
    integrateVector(WeibOSInteg,
                    upper = t,
                    h01 = h01, h02 = h02, h12 = h12, p01 = p01, p02 = p02, p12 = p12
    )
}


#' Cumulative Hazard for Piecewise Constant Hazards
#'
#' @param t (`numeric`)\cr  study time-points.
#' @param haz (`numeric vector`)\cr constant transition hazards.
#' @param pw (`numeric vector`)\cr time intervals for the piecewise constant hazards.
#'
#' @return This returns the value of cumulative hazard at time t.
#' @export
#'
#' @examples
#' pwA(1:5, c(0.5, 0.9), c(0, 4))
pwA <- function(t, haz, pw) {
  assert_numeric(t, lower = 0)
  assert_numeric(haz, lower = 0, any.missing = FALSE, all.missing = FALSE)
  assert_intervals(pw, length(haz))

  # Time intervals: time-points + cut points.
  pw_new <- c(pw[pw < max(t)])
  pw_time <- sort(unique(c(pw_new, t)))
  haz_all <- getPCWHazard(haz, pw, pw_time)

  # Cumulative hazard.
  i <- 1:(length(pw_time) - 1)
  dA_pw <- haz_all[i] * (pw_time[i + 1] - pw_time[i])
  A_pw <- c(0, cumsum(dA_pw))

  # Match result with input.
  A_pw[match(t, pw_time)]
}

#' PFS Survival Function from Piecewise Constant Hazards
#'
#' @inheritParams PWCsurvOS
#'
#' @return This returns the value of PFS survival function at time t.
#' @export
#'
#' @examples
#' PWCsurvPFS(1:5, c(0.3, 0.5), c(0.5, 0.8), c(0, 4), c(0, 8))
PWCsurvPFS <- function(t, h01, h02, pw01, pw02) {
  assert_numeric(t, lower = 0)
  assert_numeric(h01, lower = 0, any.missing = FALSE, all.missing = FALSE)
  assert_numeric(h02, lower = 0, any.missing = FALSE, all.missing = FALSE)

  assert_intervals(pw01, length(h01))
  assert_intervals(pw02, length(h02))

  A01 <- pwA(t, h01, pw01)
  A02 <- pwA(t, h02, pw02)

  exp(-A01 - A02)
}

#' Helper Function of `PWCsurvOS()`
#'
#' @param x (`numeric`)\cr  variable of integration.
#' @inheritParams PWCsurvOS
#'
#' @return Numeric results of the integrand used to calculate
#' the OS survival function for piecewise constant transition hazards, see  `PWCsurvOS`.
#' @export
#'
#' @examples
#' PwcOSInt(1:5, c(0.3, 0.5), c(0.5, 0.8), c(0.7, 1), c(0, 4), c(0, 8), c(0, 3))
PwcOSInt <- function(x, h01, h02, h12, pw01, pw02, pw12) {
  PWCsurvPFS(x, h01, h02, pw01, pw02) *
    getPCWHazard(h01, pw01, x) *
    exp(pwA(x, h12, pw12))
}

#' OS Survival Function from Piecewise Constant Hazards
#'
#' @param t (`numeric`)\cr study time-points.
#' @param h01 (`numeric vector`)\cr constant transition hazards for 0 to 1 transition.
#' @param h02 (`numeric vector`)\cr constant transition hazards for 0 to 2 transition.
#' @param h12 (`numeric vector`)\cr constant transition hazards for 1 to 2 transition.
#' @param pw01 (`numeric vector`)\cr time intervals for the piecewise constant hazards `h01`.
#' @param pw02 (`numeric vector`)\cr time intervals for the piecewise constant hazards `h02`.
#' @param pw12 (`numeric vector`)\cr time intervals for the piecewise constant hazards `h12`.
#'
#' @return This returns the value of OS survival function at time t.
#' @export

#'
#' @examples
#' PWCsurvOS(1:5, c(0.3, 0.5), c(0.5, 0.8), c(0.7, 1), c(0, 4), c(0, 8), c(0, 3))
PWCsurvOS <- function(t, h01, h02, h12, pw01, pw02, pw12) {
  assert_numeric(h01, lower = 0, any.missing = FALSE, all.missing = FALSE)
  assert_numeric(h02, lower = 0, any.missing = FALSE, all.missing = FALSE)
  assert_numeric(h12, lower = 0, any.missing = FALSE, all.missing = FALSE)
  assert_intervals(pw01, length(h01))
  assert_intervals(pw02, length(h02))
  assert_intervals(pw12, length(h12))

  PWCsurvPFS(t, h01, h02, pw01, pw02) +
    exp(-pwA(t, h12, pw12)) *
    integrateVector(PwcOSInt,
                    upper = t,
                    h01 = h01, h02 = h02, h12 = h12, pw01 = pw01, pw02 = pw02, pw12 = pw12
    )
}
