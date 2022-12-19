
#' PFS survival function from constant transition hazards.
#'
#' @param t (positive `number`)\cr  study time-point.
#' @param h01 (positive `number`)\cr transition hazard for 0 to 1 transition.
#' @param h02 (positive `number`)\cr transition hazard for 0 to 2 transition.
#'
#' @return This returns the value of PFS survival function at time t.
#' @export
#'
#' @examples
#' ExpSurvPFS(5, 0.2, 0.4)
ExpSurvPFS <- function(t, h01, h02) {
  assert_positive_number(t, zero_ok = TRUE)
  assert_positive_number(h01, zero_ok = TRUE)
  assert_positive_number(h02, zero_ok = TRUE)

  res <- exp(-(h01 + h02) * t)
  return(res)
}


#' OS survival function from constant transition hazards.
#'
#' @param t (positive `number`)\cr  study time-point.
#' @param h01  (positive `number`)\cr transition hazard for 0 to 1 transition.
#' @param h02 (positive `number`)\cr transition hazard for 0 to 2 transition.
#' @param h12 (positive `number`)\cr transition hazard for 1 to 2 transition.
#'
#' @return This returns the value of OS survival function at time t.
#' @export
#'
#' @examples
#' ExpSurvOS(5, 0.2, 0.4, 0.1)
ExpSurvOS <- function(t, h01, h02, h12) {
  assert_positive_number(t, zero_ok = TRUE)
  assert_positive_number(h01, zero_ok = TRUE)
  assert_positive_number(h02, zero_ok = TRUE)
  assert_positive_number(h12, zero_ok = TRUE)

  h012 <- h12 - h01 - h02
  res <- ExpSurvPFS(t, h01, h02) / h012 * (h12 - h02 - h01 * exp(-h012 * t))
  return(res)
}



#' PFS survival function from Weibull transition hazards.
#'
#' @param t (positive `number`)\cr  study time-point.
#' @param h01 (positive `number`)\cr transition hazard for 0 to 1 transition.
#' @param h02 (positive `number`)\cr transition hazard for 0 to 2 transition.
#' @param p01 (positive `number`)\cr rate parameter of Weibull distribution for `h01`.
#' @param p02 (positive `number`)\cr rate parameter of Weibull distribution for `h02`.
#'
#' @return  This returns the value of PFS survival function at time t.
#' @export
#' @export
#'
#' @examples
#' WeibSurvPFS(5, 0.2, 0.5, 1.2, 0.9)
WeibSurvPFS <- function(t, h01, h02, p01, p02) {
  assert_positive_number(t, zero_ok = TRUE)
  assert_positive_number(h01, zero_ok = TRUE)
  assert_positive_number(h02, zero_ok = TRUE)
  assert_positive_number(p01)
  assert_positive_number(p02)

  res <- exp(-h01 * t^p01 - h02 * t^p02)
  return(res)
}


#' OS survival function from Weibull transition hazards.
#'
#' @param t (positive `number`)\cr  study time-point.
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
#' @examples WeibSurvOS(5, 0.2, 0.5, 2.1, 1.2, 0.9, 1)
WeibSurvOS <- function(t, h01, h02, h12, p01, p02, p12) {
  assert_positive_number(t, zero_ok = TRUE)
  assert_positive_number(h01, zero_ok = TRUE)
  assert_positive_number(h02, zero_ok = TRUE)
  assert_positive_number(h12, zero_ok = TRUE)
  assert_positive_number(p01)
  assert_positive_number(p02)
  assert_positive_number(p12)

  integrand <- function(x) {
    x^(p01 - 1) * exp(-h01 * x^p01 - h02 * x^p02 - h12 * (t^p12 - x^p12))
  }
  res <- WeibSurvPFS(t, h01, h02, p01, p02) +
    h01 * p01 * integrate(integrand, lower = 0, upper = t)$value
  return(res)
}


#' Cumulative hazard for piecewise constant hazards
#'
#' @param t (positive `number`)\cr  study time-point.
#' @param haz (`numeric vector`)\cr constant transition hazards.
#' @param pw (`numeric vector`)\cr time intervals for the piecewise constant hazards.
#'
#' @return This returns the value of cumulative hazard at time t.
#' @export
#'
#' @examples
#' pwA(5, c(0.5, 0.9), c(0, 4))
pwA <- function(t, haz, pw) {
  assert_positive_number(t, zero_ok = TRUE)
  assert_numeric(haz, lower = 0, any.missing = FALSE, all.missing = FALSE)
  assert_intervals(pw, length(haz))

  # adjust time intervals
  pw_new <- c(pw[pw < t], t)
  # cumulative hazard
  if (t != 0) {
    i <- 1:(length(pw_new) - 1)
    A_pw <- sum(haz[i] * (pw_new[i + 1] - pw_new[i]))
    res <- A_pw
  } else {
    res <- 0
  }
  return(res)
}


#' PFS survival function from piecewise constant hazards.
#'
#' @param t  (positive `number`)\cr study time-point.
#' @param h01 (`numeric vector`)\cr constant transition hazards for 0 to 1 transition
#' @param h02 (`numeric vector`)\cr constant transition hazards for 0 to 2 transition
#' @param pw01 (`numeric vector`)\cr time intervals for the piecewise constant hazards `h01`.
#' @param pw02 (`numeric vector`)\cr time intervals for the piecewise constant hazards `h02`.
#'
#' @return This returns the value of PFS survival function at time t.
#' @export
#'
#' @examples
#' PWCsurvPFS(5, c(0.3, 0.5), c(0.5, 0.8), c(0, 4), c(0, 8))
PWCsurvPFS <- function(t, h01, h02, pw01, pw02) {
  assert_positive_number(t, zero_ok = TRUE)
  assert_numeric(h01, lower = 0, any.missing = FALSE, all.missing = FALSE)
  assert_numeric(h02, lower = 0, any.missing = FALSE, all.missing = FALSE)

  assert_intervals(pw01, length(h01))
  assert_intervals(pw02, length(h02))

  A01 <- pwA(t, h01, pw01)
  A02 <- pwA(t, h02, pw02)

  res <- exp(-A01 - A02)
  return(res)
}



#' OS survival function from piecewise constant hazards.
#'
#' @param t (positive `number`)\cr study time-point.
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
#' PWCsurvOS(5, c(0.3, 0.5), c(0.5, 0.8), c(0.7, 1), c(0, 4), c(0, 8), c(0, 3))
PWCsurvOS <- function(t, h01, h02, h12, pw01, pw02, pw12) {
  integrand <- function(x) {
    sapply(x, function(u) {
      return(PWCsurvPFS(u, h01, h02, pw01, pw02) *
        getPCWHazard(h01, pw01, u) *
        exp(pwA(u, h12, pw12)))
    })
  }
  res <- PWCsurvPFS(t, h01, h02, pw01, pw02) +
    exp(-pwA(t, h12, pw12)) * integrate(integrand, lower = 0, upper = t)$value

  return(res)
}
