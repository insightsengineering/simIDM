#' PFS Survival Function from Constant Transition Hazards
#'
#' @inheritParams WeibSurvOS
#' @return This returns the value of PFS survival function at time t.
#' @export
#'
#' @examples
#' ExpSurvPFS(c(1:5), 0.2, 0.4)
ExpSurvPFS <- function(t, h01, h02) {
  assert_numeric(t, lower = 0, any.missing = FALSE)
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
  assert_numeric(t, lower = 0, any.missing = FALSE)
  assert_positive_number(h01, zero_ok = TRUE)
  assert_positive_number(h02, zero_ok = TRUE)
  assert_positive_number(h12, zero_ok = TRUE)

  h012 <- h12 - h01 - h02
  ExpSurvPFS(t, h01, h02) + h01 * h012^-1 * (ExpSurvPFS(t, h01, h02) - exp(-h12 * t))
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
  assert_numeric(t, lower = 0, any.missing = FALSE)
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
#' @keywords internal
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
#'
#' @keywords internal
WeibOSInteg <- function(x, t, h01, h02, h12, p01, p02, p12) {
  assert_numeric(x, finite = TRUE, any.missing = FALSE)
  assert_numeric(t, finite = TRUE, any.missing = FALSE)
  assert_true(test_scalar(x) || identical(length(x), length(t)) || test_scalar(t))
  assert_positive_number(h01, zero_ok = TRUE)
  assert_positive_number(h02, zero_ok = TRUE)
  assert_positive_number(h12, zero_ok = TRUE)
  assert_positive_number(p01)
  assert_positive_number(p02)
  assert_positive_number(p12)

  x^(p01 - 1) * exp(-h01 * x^p01 - h02 * x^p02 - h12 * (t^p12 - x^p12))
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
  assert_numeric(t, lower = 0, any.missing = FALSE)
  assert_positive_number(h01, zero_ok = TRUE)
  assert_positive_number(p01)

  WeibSurvPFS(t, h01, h02, p01, p02) +
    h01 * p01 *
      sapply(t, function(t) {
        integrateVector(WeibOSInteg,
          upper = t,
          t = t,
          h01 = h01,
          h02 = h02,
          h12 = h12,
          p01 = p01,
          p02 = p02,
          p12 = p12
        )
      })
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
  assert_numeric(t, lower = 0, any.missing = FALSE)
  assert_numeric(haz, lower = 0, any.missing = FALSE, all.missing = FALSE)
  assert_intervals(pw, length(haz))

  # Time intervals: time-points + cut points.
  pw_new <- c(pw[pw < max(t)])
  pw_time <- sort(unique(c(pw_new, t)))
  haz_all <- getPWCHazard(haz, pw, pw_time)

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
  assert_numeric(t, lower = 0, any.missing = FALSE)
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
#'
#' @keywords internal
PwcOSInt <- function(x, t, h01, h02, h12, pw01, pw02, pw12) {
  PWCsurvPFS(x, h01, h02, pw01, pw02) *
    getPWCHazard(h01, pw01, x) *
    exp(pwA(x, h12, pw12) - pwA(t, h12, pw12))
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
  assert_numeric(t, lower = 0, any.missing = FALSE)
  assert_numeric(h01, lower = 0, any.missing = FALSE, all.missing = FALSE)
  assert_numeric(h02, lower = 0, any.missing = FALSE, all.missing = FALSE)
  assert_numeric(h12, lower = 0, any.missing = FALSE, all.missing = FALSE)
  assert_intervals(pw01, length(h01))
  assert_intervals(pw02, length(h02))
  assert_intervals(pw12, length(h12))

  # Assemble unique time points and corresponding hazard values once.
  unique_pw_times <- sort(unique(c(pw01, pw02, pw12)))
  h01_at_times <- getPWCHazard(h01, pw01, unique_pw_times)
  h02_at_times <- getPWCHazard(h02, pw02, unique_pw_times)
  h12_at_times <- getPWCHazard(h12, pw12, unique_pw_times)

  # We work for the integral parts with sorted and unique time points.
  t_sorted <- sort(unique(t))
  t_sorted_zero <- c(0, t_sorted)
  n_t <- length(t_sorted)
  int_sums <- numeric(n_t)
  cum_haz_12 <- pwA(t_sorted, h12, pw12)
  for (i in seq_len(n_t)) {
    t_start <- t_sorted_zero[i]
    t_end <- t_sorted_zero[i + 1]
    # Determine the indices of the time intervals we are working with here.
    start_index <- findInterval(t_start, unique_pw_times)
    # We use here left.open = TRUE, so that when t_end is identical to a change point,
    # then we don't need an additional part.
    end_index <- max(findInterval(t_end, unique_pw_times, left.open = TRUE), 1L)
    index_bounds <- start_index:end_index
    # Determine corresponding time intervals.
    time_bounds <- if (length(index_bounds) == 1L) {
      rbind(c(t_start, t_end))
    } else if (length(index_bounds) == 2L) {
      rbind(
        c(t_start, unique_pw_times[end_index]),
        c(unique_pw_times[end_index], t_end)
      )
    } else {
      n_ind_bounds <- length(index_bounds)
      rbind(
        c(t_start, unique_pw_times[start_index + 1L]),
        cbind(
          unique_pw_times[index_bounds[-c(1, n_ind_bounds)]],
          unique_pw_times[index_bounds[-c(1, 2)]]
        ),
        c(unique_pw_times[end_index], t_end)
      )
    }
    for (j in seq_along(index_bounds)) {
      time_j <- time_bounds[j, 1]
      time_jp1 <- time_bounds[j, 2]
      intercept_a <- pwA(time_j, h12, pw12) - pwA(time_j, h01, pw01) - pwA(time_j, h02, pw02)
      slope_b <- h12_at_times[index_bounds[j]] - h01_at_times[index_bounds[j]] - h02_at_times[index_bounds[j]]
      h01_j <- h01_at_times[index_bounds[j]]
      int_js <- if (slope_b != 0) {
        h01_j / slope_b * exp(intercept_a) *
          (exp(slope_b * (time_jp1 - time_j) - cum_haz_12[i:n_t]) - exp(-cum_haz_12[i:n_t]))
      } else {
        h01_j * exp(intercept_a - cum_haz_12[i:n_t]) * (time_jp1 - time_j)
      }
      int_sums[i:n_t] <- int_sums[i:n_t] + int_js
    }
  }

  # Match back to all time points,
  # which might not be unique or ordered.
  int_sums <- int_sums[match(t, t_sorted)]
  result <- PWCsurvPFS(t, h01, h02, pw01, pw02) + int_sums

  # Cap at 1 in a safe way - i.e. first check.
  assert_numeric(result, finite = TRUE, any.missing = FALSE)
  above_one <- result > 1
  if (any(above_one)) {
    assert_true(all((result[above_one] - 1) < sqrt(.Machine$double.eps)))
    result[above_one] <- 1
  }
  result
}

#' Helper Function for Single Quantile for OS Survival Function
#'
#' @param q (`number`)\cr single quantile at which to compute event time.
#' @inheritParams ExpQuantOS
#'
#' @return Single time t such that the OS survival function at t equals q.
#' @keywords internal
singleExpQuantOS <- function(q, h01, h02, h12) {
  assert_number(q, finite = TRUE, lower = 0, upper = 1)
  toroot <- function(x) {
    return(q - ExpSurvOS(t = x, h01, h02, h12))
  }
  stats::uniroot(
    toroot,
    interval = c(0, 10^3),
    extendInt = "yes"
  )$root
}

#' Quantile function for OS survival function induced by an illness-death model
#'
#' @param q (`numeric`)\cr quantile(s) at which to compute event time (q = 1 / 2 corresponds to median).
#' @param h01 (`numeric vector`)\cr constant transition hazards for 0 to 1 transition.
#' @param h02 (`numeric vector`)\cr constant transition hazards for 0 to 2 transition.
#' @param h12 (`numeric vector`)\cr constant transition hazards for 1 to 2 transition.
#'
#' @return This returns the time(s) t such that the OS survival function at t equals q.
#' @export
#'
#' @examples ExpQuantOS(1 / 2, 0.2, 0.5, 2.1)
ExpQuantOS <- function(q = 1 / 2, h01, h02, h12) {
  assert_numeric(q, lower = 0, upper = 1, any.missing = FALSE)
  assert_numeric(h01, lower = 0, any.missing = FALSE)
  assert_numeric(h02, lower = 0, any.missing = FALSE)
  assert_numeric(h12, lower = 0, any.missing = FALSE)

  vapply(
    q,
    FUN = singleExpQuantOS,
    FUN.VALUE = numeric(1),
    h01 = h01,
    h02 = h02,
    h12 = h12
  )
}
