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
#' - trans (`integer`): transition (1, 2 or 3) identifier.
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
  with(data, -sum(log(haz(transition, exit, trans)^status
    * survTrans(transition, exit, trans) / survTrans(transition, entry, trans))))
}

#' Hazard Function for Different Transition Models
#'
#' @param transition (`ExponentialTransition` or `WeibullTransition`)\cr
#'   see [exponential_transition()] or [weibull_transition()] for details.
#' @param t (`numeric`)\cr time at which hazard is to be computed.
#' @param trans (`integer`)\cr index specifying the transition type.
#'
#' @return The hazard rate for the specified transition and time.
#' @export
#'
#' @details
#' This function dispatches to either `haz.ExponentialTransition` or `haz.WeibullTransition`
#' based on the `transition` object class.
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
  params[trans]
}

#' Hazard Function for Weibull Transition Model
#'
#' @param transition (`WeibullTransition`)\cr
#'   See [weibull_transition()] for details.
#' @param t (`numeric`)\cr time at which hazard is to be computed.
#' @param trans (`integer`)\cr index specifying the transition type.
#'
#' @return Returns the hazard rate (`numeric`) for the Weibull transition at the specified time.
#' @export
#'
#' @details
#' Computes the hazard function specifically for a Weibull transition model using provided parameters.
#'
#' @examples
#' transition <- weibull_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6, p01 = 2, p02 = 2.5, p12 = 3)
#' haz(transition, 0.4, 2)
haz.WeibullTransition <- function(transition, t, trans) {
  # params (in this order): h01, h02, h12, p01, p02, p12.
  params <- c(unlist(transition$hazards, use.names = FALSE), unlist(transition$weibull_rates, use.names = FALSE))
  params[trans] * params[trans + 3] * t^(params[trans + 3] - 1)
}

#' Survival Function for Different Transition Models
#'
#' @param transition (`ExponentialTransition` or `WeibullTransition`)\cr
#'   See [exponential_transition()] or [weibull_transition()] for details.
#' @param t (`numeric`)\cr time at which survival probability is to be computed.
#' @param trans (`integer`)\cr index specifying the transition type.
#'
#' @return Returns the survival probability (`numeric`) for the specified transition and time.
#' @export
#'
#' @details
#' This function dispatches to either `survTrans.ExponentialTransition` or `survTrans.WeibullTransition`
#' based on the `transition` object class
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

#' Survival Function for Exponential Transition Model
#'
#' @param transition (`ExponentialTransition`)\cr
#'   See [exponential_transition()] for details.
#' @param t (`numeric`)\cr time at which survival probability is to be computed.
#' @param trans (`integer`)\cr index specifying the transition type.
#'
#' @return Returns the survival probability (`numeric`) for the exponential transition at the specified time.
#' @export
#'
#' @details
#' Computes the survival function specifically for an exponential transition model using provided parameters.
#'
#' @examples
#' transition <- exponential_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6)
#' survTrans(transition, 0.4, 2)
survTrans.ExponentialTransition <- function(transition, t, trans) {
  # params (in this order): h01, h02, h12.
  params <- unlist(transition$hazards, use.names = FALSE)
  exp(-params[trans] * t)
}

#' Survival Function for Weibull Transition Model
#'
#' @param transition (`WeibullTransition`)\cr
#'   See [weibull_transition()] for details.
#' @param t (`numeric`)\cr time at which survival probability is to be computed.
#' @param trans (`integer`)\cr index specifying the transition type.
#'
#' @return Returns the survival probability (`numeric`) for the Weibull transition at the specified time.
#' @export
#'
#' @details
#' Computes the survival function specifically for a Weibull transition model using provided parameters.
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
#' @return Returns a vector (`numeric`) of initial parameters for likelihood maximization.
#' @export
#'
#' @details
#' This function dispatches to either `getInit.ExponentialTransition` or `getInit.WeibullTransition`
#' based on the `transition` object class.
#'
#' @examples
#' transition <- exponential_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6)
#' getInit(transition)
getInit <- function(transition) {
  assert_class(transition, "TransitionParameters")
  UseMethod("getInit")
}

#' Retrieve Initial Parameters for Exponential Transition Model
#'
#' @param transition (`ExponentialTransition`)\cr containing the initial parameters.
#'   See [exponential_transition()] for details.
#'
#' @return Returns a vector (`numeric`) of initial parameters for the exponential transition.
#' @export
#'
#' @details
#' Extracts initial parameter values specifically for an exponential transition model.
#'
#' @examples
#' transition <- exponential_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6)
#' getInit(transition)
getInit.ExponentialTransition <- function(transition) {
  unlist(transition$hazards, use.names = FALSE)
}

#' Retrieve Initial Parameters for Weibull Transition Model
#'
#' @param transition (`WeibullTransition`)\cr containing the initial parameters.
#'   See [weibull_transition()] for details.
#'
#' @return Returns a vector (`numeric`) of initial parameters for the Weibull transition.
#' @export
#'
#' @details
#' Extracts initial parameter values specifically for a Weibull transition model.
#'
#' @examples
#' transition <- weibull_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6, p01 = 2, p02 = 2.5, p12 = 3)
#' getInit(transition)
getInit.WeibullTransition <- function(transition) {
  c(unlist(transition$hazards, use.names = FALSE), unlist(transition$weibull_rates, use.names = FALSE))
}

#' Generate the Target Function for Optimization
#'
#' @param params (`numeric` vector)\cr
#'   Vector of parameters to be optimized over.
#' @param data (`data.frame`)\cr
#'   Data frame in the format created by [prepareData()].
#' @param transition (`transitionParameters` object)\cr
#'   Object specifying the distribution family. See [exponential_transition()] or [weibull_transition()] for details.
#'
#' @return Returns a function that calculates the negative log-likelihood (`numeric`) for the given parameters.
#' @export
#'
#' @details
#' This function creates a target function for optimization, which computes the negative log-likelihood for given
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
#' getTarget(params, data, transition)
getTarget <- function(params, data, transition) {
  assert_numeric(params)
  assert_data_frame(data)
  assert_class(transition, "TransitionParameters")
  if ("ExponentialTransition" %in% class(transition)) {
    negLogLik(transition = exponential_transition(h01 = params[1], h02 = params[2], h12 = params[3]), data = data)
  } else {
    negLogLik(transition = weibull_transition(
      h01 = params[1], h02 = params[2], h12 = params[3],
      p01 = params[4], p02 = params[5], p12 = params[6]
    ), data = data)
  }
}

#' Format Results of Parameter Estimation for Different Transition Models
#'
#' @param transition (`ExponentialTransition` or `WeibullTransition`)\cr
#'   See [exponential_transition()] or [weibull_transition()] for details.
#' @param res (`numeric` vector)\cr vector of parameter estimates from the likelihood maximization procedure.
#'
#' @return Returns a `TransitionParameters` object with parameter estimates.
#' @export
#'
#' @details
#' This function dispatches to either `getResults.ExponentialTransition` or `getResults.WeibullTransition`
#' based on the `transition` object class, using the results from likelihood maximization.
#'
#' @examples
#' results <- c(1.2, 1.5, 1.6)
#' getResults(exponential_transition(), results)
getResults <- function(transition, res) {
  UseMethod("getResults")
}

#' Format Results of Parameter Estimation for Exponential Transition
#'
#' @param transition (`ExponentialTransition`)\cr
#'   See [exponential_transition()] for details.
#' @inheritParams getResults
#'
#' @return Returns an `ExponentialTransition` object with updated parameter estimates.
#' @export
#'
#' @details
#' Constructs an `ExponentialTransition` object with parameters estimated from likelihood maximization.
#'
#' @examples
#' results <- c(1.2, 1.5, 1.6)
#' getResults(exponential_transition(), results)
getResults.ExponentialTransition <- function(transition, res) {
  exponential_transition(h01 = res[1], h02 = res[2], h12 = res[3])
}

#' Format Results of Parameter Estimation for Weibull Transition
#'
#' @param transition (`WeibullTransition`)\cr
#'   See [weibull_transition()] for details.
#' @inheritParams getResults
#'
#' @return Returns a `WeibullTransition` object with updated parameter estimates.
#' @export
#'
#' @details
#' Constructs a `WeibullTransition` object with parameters estimated from likelihood maximization.
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
#' @param data (`data.frame`)\cr
#'   Data frame in the format produced by [getOneClinicalTrial()].
#' @param transition (`transitionParameters` object)\cr
#'   Object specifying the assumed distribution of transition hazards.
#'   Initial parameters for optimization can be specified here.
#'   See [exponential_transition()] or [weibull_transition()] for details.
#'
#' @return Returns a `transitionParameters` object with the estimated parameters.
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

  res <- stats::optim(
    par = getInit(transition),
    fn = getTarget,
    method = "L-BFGS-B",
    lower = 1e-3,
    data = data,
    transition = transition
  )$par

  getResults(transition, res)
}
