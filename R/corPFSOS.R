# survPFS ----

#' PFS Survival Function for Different Transition Models
#'
#' @param transition (`TransitionParameters`)\cr
#'   see [exponential_transition()], [weibull_transition()] or [piecewise_exponential()] for details.
#' @param t (`numeric`)\cr time at which the value of the PFS survival function is to be computed.
#'
#' @return The value of the survival function for the specified transition and time.
#' @export
#'
#' @examples
#' transition <- exponential_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6)
#' survPFS(transition, 0.4)
survPFS <- function(transition, t) {
  UseMethod("survPFS")
}

#' @describeIn survPFS Survival Function for an exponential transition model.
#' @export
#'
#' @examples
#' transition <- exponential_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6)
#' survPFS(transition, 0.4)
survPFS.ExponentialTransition <- function(transition, t) {
  ExpSurvPFS(t = t, h01 = transition$hazards$h01, h02 = transition$hazards$h02)
}

#' @describeIn survPFS Survival Function for a Weibull transition model.
#' @export
#'
#' @examples
#' transition <- weibull_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6, p01 = 2, p02 = 2.5, p12 = 3)
#' survPFS(transition, 0.4)
survPFS.WeibullTransition <- function(transition, t) {
  WeibSurvPFS(
    t = t, h01 = transition$hazards$h01, h02 = transition$hazards$h02,
    p01 = transition$weibull_rates$p01, p02 = transition$weibull_rates$p02
  )
}

#' @describeIn survPFS Survival Function for a piecewise constant transition model.
#' @export
#'
#' @examples
#' transition <- piecewise_exponential(
#'   h01 = c(1, 1, 1), h02 = c(1.5, 0.5, 1), h12 = c(1, 1, 1),
#'   pw01 = c(0, 3, 8), pw02 = c(0, 6, 7), pw12 = c(0, 8, 9)
#' )
#' survPFS(transition, 0.4)
survPFS.PWCTransition <- function(transition, t) {
  PWCsurvPFS(t,
    h01 = transition$hazards$h01, h02 = transition$hazards$h02,
    pw01 = transition$intervals$pw01, pw02 = transition$intervals$pw02
  )
}

# survOS ----

#' OS Survival Function for Different Transition Models
#'
#' @param transition (`TransitionParameters`)\cr
#'   see [exponential_transition()], [weibull_transition()] or [piecewise_exponential()] for details.
#' @param t (`numeric`)\cr time at which the value of the OS survival function is to be computed.
#'
#' @return The value of the survival function for the specified transition and time.
#' @export
#'
#' @examples
#' transition <- exponential_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6)
#' survOS(transition, 0.4)
survOS <- function(transition, t) {
  UseMethod("survOS")
}

#' @describeIn survOS Survival Function for an exponential transition model.
#' @export
#'
#' @examples
#' transition <- exponential_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6)
#' survOS(transition, 0.4)
survOS.ExponentialTransition <- function(transition, t) {
  ExpSurvOS(
    t = t, h01 = transition$hazards$h01, h02 = transition$hazards$h02,
    h12 = transition$hazards$h12
  )
}

#' @describeIn survOS Survival Function for a Weibull transition model.
#' @export
#'
#' @examples
#' transition <- weibull_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6, p01 = 2, p02 = 2.5, p12 = 3)
#' survOS(transition, 0.4)
survOS.WeibullTransition <- function(transition, t) {
  WeibSurvOS(
    t = t, h01 = transition$hazards$h01, h02 = transition$hazards$h02,
    h12 = transition$hazards$h12, p01 = transition$weibull_rates$p01,
    p02 = transition$weibull_rates$p02, p12 = transition$weibull_rates$p12
  )
}

#' @describeIn survOS Survival Function for a piecewise constant transition model.
#' @export
#'
#' @examples
#' transition <- piecewise_exponential(
#'   h01 = c(1, 1, 1), h02 = c(1.5, 0.5, 1), h12 = c(1, 1, 1),
#'   pw01 = c(0, 3, 8), pw02 = c(0, 6, 7), pw12 = c(0, 8, 9)
#' )
#' survOS(transition, 0.4)
survOS.PWCTransition <- function(transition, t) {
  PWCsurvOS(
    t = t, h01 = transition$hazards$h01, h02 = transition$hazards$h02,
    h12 = transition$hazards$h12, pw01 = transition$intervals$pw01,
    pw02 = transition$intervals$pw02, pw12 = transition$intervals$pw12
  )
}

# expval ----

#' Helper Function for Computing E(PFS^2)
#'
#' @param x (`numeric`)\cr variable of integration.
#' @param transition (`TransitionParameters`)\cr
#'   see [exponential_transition()], [weibull_transition()] or [piecewise_exponential()] for details.
#'
#' @return Numeric results of the integrand used to calculate E(PFS^2).
#' @export
#'
#' @examples
#' transition <- exponential_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6)
#' expvalPFSInteg(0.4, transition)
expvalPFSInteg <- function(x, transition) {
  x * survPFS(transition, x)
}

#' Helper Function for Computing E(OS^2)
#'
#' @param x (`numeric`)\cr variable of integration.
#' @param transition (`TransitionParameters`)\cr
#'   see [exponential_transition()], [weibull_transition()] or [piecewise_exponential()] for details.
#'
#' @return Numeric results of the integrand used to calculate E(OS^2).
#' @export
#'
#' @examples
#' transition <- exponential_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6)
#' expvalOSInteg(0.4, transition)
expvalOSInteg <- function(x, transition) {
  x * survOS(transition = transition, t = x)
}

# p11 ----

#' Helper Function for `log_p11()`
#'
#' @param x (`numeric`)\cr variable of integration.
#' @param transition (`TransitionParameters`)\cr
#'   see [exponential_transition()], [weibull_transition()] or [piecewise_exponential()] for details.
#'
#' @return Hazard rate at the specified time for the transition from progression to death.
#' @keywords internal
p11Integ <- function(x, transition) {
  haz(transition = transition, t = x, trans = 3)
}

#' Probability of Remaining in Progression Between Two Time Points for Different Transition Models
#'
#' @param transition (`TransitionParameters`)\cr
#'   see [exponential_transition()], [weibull_transition()] or [piecewise_exponential()] for details.
#' @param s (`numeric`)\cr lower time points.
#' @param t (`numeric`)\cr higher time points.
#' @return This returns the natural logarithm of the probability of remaining in progression (state 1)
#' between two time points, conditional on being in state 1 at the lower time point.
#'
#' @export
#'
#' @examples
#' transition <- exponential_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6)
#' log_p11(transition, 1, 3)
log_p11 <- function(transition, s, t) {
  assert_numeric(s, finite = TRUE, any.missing = FALSE, lower = 0)
  assert_numeric(t, finite = TRUE, any.missing = FALSE, lower = 0)
  assert_true(identical(length(s), length(t)))
  assert_true(all(t > s))

  intval <- mapply(function(s, t) {
    stats::integrate(p11Integ,
      lower = s,
      upper = t,
      transition
    )$value
  }, s, t)
  -intval
}

# PFSOS ----

#' Helper Function for `survPFSOS()`
#'
#' @param u (`numeric`)\cr variable of integration.
#' @param t (`numeric`)\cr time at which the value of the PFS*OS survival function is to be computed.
#' @param transition (`TransitionParameters`)\cr
#'   see [exponential_transition()], [weibull_transition()] or [piecewise_exponential()] for details.
#'
#' @note Not all vectors `u` and `t` work here due to assertions in [log_p11()].
#'
#' @return Numeric result of the integrand used to calculate the PFS*OS survival function.
#' @keywords internal
PFSOSInteg <- function(u, t, transition) {
  exp(log_p11(transition, u, t / u) + log(survPFS(transition, u)) + log(haz(transition, u, 1)))
}

#' Survival Function of the Product PFS*OS for Different Transition Models
#'
#' @param t (`numeric`)\cr time at which the value of the PFS*OS survival function is to be computed.
#' @param transition (`TransitionParameters`)\cr
#'   see [exponential_transition()], [weibull_transition()] or [piecewise_exponential()] for details.
#'
#' @return This returns the value of PFS*OS survival function at time t.
#' @export
#'
#' @examples
#' transition <- exponential_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6)
#' survPFSOS(0.4, transition)
survPFSOS <- function(t, transition) {
  sapply(t, function(x) {
    intval <- stats::integrate(PFSOSInteg, lower = 0, upper = sqrt(x), x, transition)$value
    survPFS(transition, sqrt(x)) + intval
  })
}

# correlation ----

#' Correlation of PFS and OS event times for Different Transition Models
#'
#' @param transition (`TransitionParameters`)\cr
#'   see [exponential_transition()], [weibull_transition()] or [piecewise_exponential()] for details.
#'
#' @return The correlation of PFS and OS.
#' @export
#'
#' @examples
#' transition <- exponential_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6)
#' corTrans(transition)
corTrans <- function(transition) {
  # E(PFS) & E(OS).
  expvalPFS <- stats::integrate(survPFS,
    lower = 0, upper = Inf,
    transition = transition
  )$value

  expvalOS <- stats::integrate(survOS,
    lower = 0, upper = Inf,
    transition = transition
  )$value

  # Var(PFS) & Var(OS).
  expvalPFS2 <- 2 * stats::integrate(expvalPFSInteg,
    lower = 0, upper = Inf,
    transition = transition
  )$value

  expvalOS2 <- 2 * stats::integrate(expvalOSInteg,
    lower = 0, upper = Inf,
    transition = transition
  )$value

  varPFS <- expvalPFS2 - expvalPFS^2

  varOS <- expvalOS2 - expvalOS^2

  # E(PFS*OS).
  expvalPFSOS <- stats::integrate(survPFSOS,
    lower = 0, upper = Inf,
    transition
  )$value

  # Cor(PFS, OS).
  (expvalPFSOS - expvalPFS * expvalOS) / sqrt(varPFS * varOS)
}

#' Correlation of PFS and OS event times for data from the IDM
#'
#' @param data (`data.frame`)\cr in the format produced by [getOneClinicalTrial()].
#' @param transition (`TransitionParameters` object)\cr specifying the assumed distribution of transition hazards.
#'   Initial parameters for optimization can be specified here.
#'   See [exponential_transition()] or [weibull_transition()] for details.
#' @param bootstrap (`flag`)\cr if `TRUE` computes confidence interval via bootstrap.
#' @param bootstrap_n (`count`)\cr number of bootstrap samples.
#' @param conf_level (`proportion`)\cr confidence level for the confidence interval.
#'
#' @return The correlation of PFS and OS.
#' @export
#'
#' @examples
#' transition <- exponential_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6)
#' data <- getClinicalTrials(
#'   nRep = 1, nPat = c(100), seed = 1234, datType = "1rowTransition",
#'   transitionByArm = list(transition), dropout = list(rate = 0.5, time = 12),
#'   accrual = list(param = "intensity", value = 7)
#' )[[1]]
#' corPFSOS(data, transition = exponential_transition(), bootstrap = FALSE)
#' \dontrun{
#' corPFSOS(data, transition = exponential_transition(), bootstrap = TRUE)
#' }
corPFSOS <- function(data, transition, bootstrap = TRUE, bootstrap_n = 100, conf_level = 0.95) {
  assert_data_frame(data)
  assert_flag(bootstrap)
  assert_count(bootstrap_n)
  assert_number(conf_level, lower = 0.01, upper = 0.999)

  trans <- estimateParams(data, transition)
  res <- list("corPFSOS" = corTrans(trans))
  if (bootstrap) {
    future::plan(future::multisession, workers = max(1, parallelly::availableCores() - 1))
    ids <- lapply(1:bootstrap_n, function(x) sample(seq_len(nrow(data)), nrow(data), replace = TRUE))
    corBootstrap <- furrr::future_map_dbl(ids, ~ {
      furrr::furrr_options(
        globals = list(data = data, transition = transition),
        packages = c("simIDM")
      )
      b_sample <- data[.x, , drop = FALSE]
      b_transition <- estimateParams(b_sample, transition)
      corTrans(b_transition)
    })
    lowerQuantile <- (1 - conf_level) / 2
    upperQuantile <- lowerQuantile + conf_level
    c(stats::quantile(corBootstrap, lowerQuantile),
      "corPFSOS" = res,
      stats::quantile(corBootstrap, upperQuantile)
    )
    res$lower <- stats::quantile(corBootstrap, lowerQuantile)
    res$upper <- stats::quantile(corBootstrap, upperQuantile)
  }
  res
}
