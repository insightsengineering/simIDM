#' Preparation of a Data Set to Compute Log-likelihood
#'
#' @param data (`data.frame`)\cr containing entry and exit times of an illness-death model.
#'   See [getOneClinicalTrial()] for details.
#'
#' @return This function returns a data set with one row per patient and transition, when the patient is at risk.
#' @export
#'
#' @details
#' The output data set contains the following columns:
#' - id (`integer`): patient id.
#' - from (`integer`): start event state.
#' - to (`integer`): end event state.
#' - trans (`integer`): transition (1, 2 or 3) identifier
#'   - `1`: Transition from state 0 (stable) to 1 (progression).
#'   - `2`: Transition from state 0 (stable) to 2 (death).
#'   - `3`: Transition from state 1 (progression) to 2 (death).
#' - entry (`numeric`): time at which the patient begins to be at risk for the transition.
#' - exit (`numeric`): time at which the patient ends to be at risk for the transition.
#' - status (`logical`): event indicator for the transition.
#'
#' @examples
#' transition <- exponential_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6)
#' simData <- getOneClinicalTrial(
#'   nPat = c(30), transitionByArm = list(transition),
#'   dropout = list(rate = 0.8, time = 12),
#'   accrual = list(param = "time", value = 1)
#' )
#' prepareData(simData)
prepareData <- function(data) {
  assert_data_frame(data, min.cols = 9, max.cols = 11)
  colNames <- c(
    "id", "trt", "PFStime", "CensoredPFS", "PFSevent", "OStime",
    "CensoredOS", "OSevent", "recruitTime",
    "OStimeCal", "PFStimeCal"
  )
  if (!all(names(data) %in% colNames)) {
    data <- getDatasetWideFormat(data)
  }

  # Transform simIDM trial data to log-likelihood-compatible format.
  # Suppress warning about how msprep handles 1 -> 3 transitions.
  dataNew <- suppressWarnings(mstate::msprep(
    time = c("recruitTime", "PFStime", "OStime"),
    status = c("trt", "PFSevent", "OSevent"),
    data = data,
    trans = mstate::trans.illdeath(),
    id = data$id
  ))
  cols <- which(names(dataNew) %in% c("Tstart", "Tstop"))
  names(dataNew)[cols] <- c("entry", "exit")
  # Correct msprep results for uncensored PFS=OS events.
  ids <- data$id[data$PFStimeCal == data$OStimeCal & data$CensoredPFS == 0]
  dataNew <- dataNew[!(dataNew$id %in% ids & dataNew$trans == 3), ]
  dataNew$status[dataNew$id %in% ids] <- abs(dataNew$status[dataNew$id %in% ids] - 1)

  as.data.frame(dataNew[, -which(names(dataNew) == "time")], row.names = seq_len(nrow(dataNew)))
}

#' Compute the Negative Log-Likelihood for a Given Data Set and Transition Model
#'
#' @param transition (`ExponentialTransition` or `WeibullTransition`)\cr
#'   see [exponential_transition()] or [weibull_transition()] for details.
#' @param data (`data.frame`)\cr in the format created by [prepareData()].
#'
#' @return The value of the negative log-likelihood.
#' @export
#'
#' @details
#' Calculates the negative log-likelihood for a given data set and transition model. It uses the hazard
#' and survival functions specific to the transition model.
#'
#' @examples
#' transition <- exponential_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6)
#' simData <- getOneClinicalTrial(
#'   nPat = c(30), transitionByArm = list(transition),
#'   dropout = list(rate = 0.8, time = 12),
#'   accrual = list(param = "time", value = 1)
#' )
#' negLogLik(transition, prepareData(simData))
negLogLik <- function(transition, data) {
  with(
    data,
    -sum(log(
      haz(transition, exit, trans)^status * survTrans(transition, exit, trans) /
        survTrans(transition, entry, trans)
    ))
  )
}

#' Hazard Function for Different Transition Models
#'
#' @param transition (`ExponentialTransition` or `WeibullTransition`)\cr
#'   see [exponential_transition()] or [weibull_transition()] for details.
#' @param t (`numeric`)\cr time at which hazard is to be computed.
#' @param trans (`integer`)\cr index specifying the transition type.
#'
#' @details
#' The transition types are:
#' - `1`: Transition from state 0 (stable) to 1 (progression).
#' - `2`: Transition from state 0 (stable) to 2 (death).
#' - `3`: Transition from state 1 (progression) to 2 (death).
#'
#' @return The hazard rate for the specified transition and time.
#' @export
#'
#' @examples
#' transition <- exponential_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6)
#' haz(transition, 0.4, 2)
haz <- function(transition, t, trans) {
  assert_class(transition, "TransitionParameters")
  assert_numeric(t, lower = 0)
  assert_subset(trans, c(1, 2, 3))
  UseMethod("haz")
}

#' @describeIn haz for an exponential transition model.
#' @export
#'
#' @examples
#' transition <- exponential_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6)
#' haz(transition, 0.4, 2)
haz.ExponentialTransition <- function(transition, t, trans) {
  # params (in this order): h01, h02, h12.
  params <- unlist(transition$hazards, use.names = FALSE)
  rep(params[trans], length(t) / length(trans))
}

#' @describeIn haz for the Weibull transition model.
#' @export
#'
#' @examples
#' transition <- weibull_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6, p01 = 2, p02 = 2.5, p12 = 3)
#' haz(transition, 0.4, 2)
haz.WeibullTransition <- function(transition, t, trans) {
  # params (in this order): h01, h02, h12, p01, p02, p12.
  params <- c(unlist(transition$hazards, use.names = FALSE), unlist(transition$weibull_rates, use.names = FALSE))
  params[trans] * params[trans + 3] * t^(params[trans + 3] - 1)
}

#' @describeIn haz for the piecewise constant transition model.
#' @export
#'
#' @examples
#' transition <- piecewise_exponential(
#'   h01 = c(1, 1, 1), h02 = c(1.5, 0.5, 1), h12 = c(1, 1, 1),
#'   pw01 = c(0, 3, 8), pw02 = c(0, 6, 7), pw12 = c(0, 8, 9)
#' )
#' haz(transition, 6, 2)
haz.PWCTransition <- function(transition, t, trans) {
  getPWCHazard(unlist(transition$hazards[trans], use.names = FALSE),
    unlist(transition$intervals[trans], use.names = FALSE),
    x = t
  )
}

#' Survival Function for Different Transition Models
#'
#' @param transition (`ExponentialTransition` or `WeibullTransition`)\cr
#'   see [exponential_transition()] or [weibull_transition()] for details.
#' @param t (`numeric`)\cr time at which survival probability is to be computed.
#' @param trans (`integer`)\cr index specifying the transition type.
#'
#' @return The survival probability for the specified transition and time.
#' @export
#'
#' @examples
#' transition <- exponential_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6)
#' survTrans(transition, 0.4, 2)
survTrans <- function(transition, t, trans) {
  assert_class(transition, "TransitionParameters")
  assert_numeric(t, lower = 0)
  assert_subset(trans, c(1, 2, 3))
  UseMethod("survTrans")
}

#' @describeIn survTrans for the Exponential Transition Model
#' @export
#'
#' @examples
#' transition <- exponential_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6)
#' survTrans(transition, 0.4, 2)
survTrans.ExponentialTransition <- function(transition, t, trans) {
  # params (in this order): h01, h02, h12.
  params <- unlist(transition$hazards, use.names = FALSE)
  exp(-params[trans] * t)
}

#' @describeIn survTrans for the Weibull Transition Model
#' @export
#'
#' @examples
#' transition <- weibull_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6, p01 = 2, p02 = 2.5, p12 = 3)
#' survTrans(transition, 0.4, 2)
survTrans.WeibullTransition <- function(transition, t, trans) {
  # params (in this order): h01, h02, h12, p01, p02, p12.
  params <- c(unlist(transition$hazards, use.names = FALSE), unlist(transition$weibull_rates, use.names = FALSE))
  exp(-params[trans] * t^params[trans + 3])
}

#' Retrieve Initial Parameter Vectors for Likelihood Maximization
#'
#' @param transition (`ExponentialTransition` or `WeibullTransition`)\cr containing the initial parameters.
#'   See [exponential_transition()] or [weibull_transition()] for details.
#'
#' @return The numeric vector of initial parameters for likelihood maximization.
#' @export
#'
#' @examples
#' transition <- exponential_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6)
#' getInit(transition)
getInit <- function(transition) {
  assert_class(transition, "TransitionParameters")
  UseMethod("getInit")
}

#' @describeIn getInit for the Exponential Transition Model
#' @export
#'
#' @examples
#' transition <- exponential_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6)
#' getInit(transition)
getInit.ExponentialTransition <- function(transition) {
  unlist(transition$hazards, use.names = FALSE)
}

#' @describeIn getInit for the Weibull Transition Model
#' @export
#'
#' @examples
#' transition <- weibull_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6, p01 = 2, p02 = 2.5, p12 = 3)
#' getInit(transition)
getInit.WeibullTransition <- function(transition) {
  c(unlist(transition$hazards, use.names = FALSE), unlist(transition$weibull_rates, use.names = FALSE))
}

#' Generate the Target Function for Optimization
#'
#' @param transition (`TransitionParameters`)\cr
#'   specifying the distribution family. See [exponential_transition()] or [weibull_transition()] for details.
#'
#' @return Function that calculates the negative log-likelihood for the given parameters.
#' @export
#'
#' @details
#' This function creates a target function for optimization, computing the negative log-likelihood for given
#' parameters, data, and transition model type.
#'
#' @examples
#' transition <- exponential_transition(2, 1.3, 0.8)
#' simData <- getOneClinicalTrial(
#'   nPat = c(30), transitionByArm = list(transition),
#'   dropout = list(rate = 0.8, time = 12),
#'   accrual = list(param = "time", value = 1)
#' )
#' params <- c(1.2, 1.5, 1.6) # For ExponentialTransition
#' data <- prepareData(simData)
#' transition <- exponential_transition()
#' fun <- getTarget(transition)
#' fun(params, data)
getTarget <- function(transition) {
  assert_class(transition, "TransitionParameters")
  UseMethod("getTarget", transition)
}

#' @describeIn getTarget for the Exponential Transition Model
#' @export
#'
#' @examples
#' transition <- exponential_transition(2, 1.3, 0.8)
#' simData <- getOneClinicalTrial(
#'   nPat = c(30), transitionByArm = list(transition),
#'   dropout = list(rate = 0.8, time = 12),
#'   accrual = list(param = "time", value = 1)
#' )
#' params <- c(1.2, 1.5, 1.6)
#' data <- prepareData(simData)
#' transition <- exponential_transition()
#' target <- getTarget(transition)
#' target(params, data)
getTarget.ExponentialTransition <- function(transition) {
  function(params, data) {
    negLogLik(transition = exponential_transition(h01 = params[1], h02 = params[2], h12 = params[3]), data = data)
  }
}

#' @describeIn getTarget for the Weibull Transition Model
#' @export
#'
#' @examples
#' transition <- weibull_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6, p01 = 2, p02 = 2.5, p12 = 3)
#' simData <- getOneClinicalTrial(
#'   nPat = c(30), transitionByArm = list(transition),
#'   dropout = list(rate = 0.8, time = 12),
#'   accrual = list(param = "time", value = 1)
#' )
#' params <- c(1.2, 1.5, 1.6, 0.8, 1.3, 1.1)
#' data <- prepareData(simData)
#' transition <- weibull_transition()
#' target <- getTarget(transition)
#' target(params, data)
getTarget.WeibullTransition <- function(transition) {
  function(params, data) {
    negLogLik(transition = weibull_transition(
      h01 = params[1], h02 = params[2], h12 = params[3],
      p01 = params[4], p02 = params[5], p12 = params[6]
    ), data = data)
  }
}

#' Format Results of Parameter Estimation for Different Transition Models
#'
#' @param transition (`TransitionParameters`)\cr
#'   see [exponential_transition()] or [weibull_transition()] for details.
#' @param res (`numeric` vector)\cr vector of parameter estimates from the likelihood maximization procedure.
#'
#' @return Returns a `TransitionParameters` object with parameter estimates.
#' @export
#'
#' @examples
#' results <- c(1.2, 1.5, 1.6)
#' getResults(exponential_transition(), results)
getResults <- function(transition, res) {
  UseMethod("getResults")
}

#' @describeIn getResults for the Exponential Transition Model
#' @export
#'
#' @examples
#' results <- c(1.2, 1.5, 1.6)
#' getResults(exponential_transition(), results)
getResults.ExponentialTransition <- function(transition, res) {
  exponential_transition(h01 = res[1], h02 = res[2], h12 = res[3])
}

#' @describeIn getResults for the Weibull Transition Model
#' @export
#'
#' @examples
#' results <- c(1.2, 1.5, 1.6, 2, 2.5, 1)
#' getResults(weibull_transition(), results)
getResults.WeibullTransition <- function(transition, res) {
  weibull_transition(
    h01 = res[1], h02 = res[2], h12 = res[3],
    p01 = res[4], p02 = res[5], p12 = res[6]
  )
}

#' Estimate Parameters of the Multistate Model Using Clinical Trial Data
#'
#' @param data (`data.frame`)\cr in the format produced by [getOneClinicalTrial()].
#' @param transition (`TransitionParameters` object)\cr specifying the assumed distribution of transition hazards.
#'   Initial parameters for optimization can be specified here.
#'   See [exponential_transition()] or [weibull_transition()] for details.
#'
#' @return Returns a `TransitionParameters` object with the estimated parameters.
#' @export
#'
#' @details
#' This function estimates parameters for transition models using clinical trial data.
#' The `transition` object can be initialized with starting values for parameter estimation.
#' It uses [stats::optim()] to optimize the parameters.
#'
#' @examples
#' transition <- exponential_transition(h01 = 2, h02 = 1.4, h12 = 1.6)
#' simData <- getOneClinicalTrial(
#'   nPat = c(30), transitionByArm = list(transition),
#'   dropout = list(rate = 0.3, time = 12),
#'   accrual = list(param = "time", value = 1)
#' )
#' # Initialize transition with desired starting values for optimization:
#' transition <- exponential_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6)
#' estimate <- estimateParams(simData, transition)
estimateParams <- function(data, transition) {
  data <- prepareData(data)
  par <- getInit(transition)
  target <- getTarget(transition)

  res <- stats::optim(
    par = par,
    fn = target,
    method = "L-BFGS-B",
    lower = 1e-3,
    data = data
  )$par

  getResults(transition, res)
}
