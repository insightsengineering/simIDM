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

negLogLik <- function(transition, data) {
  with(data, -sum(log(haz(transition, exit, trans)^status * survTrans(transition, exit, trans) / survTrans(transition, entry, trans))))
}

#' Hazard Function for Different Transition Models
#'
#' @param transition (`ExponentialTransition` or `WeibullTransition`)\cr.
#'   See [exponential_transition()] or [weibull_transition()] for details.
#' @param t (`numeric`)\cr time at which hazard is to be computed.
#' @param trans (`integer`)\cr index specifying the transition type.
#'
#' @return Returns the hazard rate (`numeric`) for the specified transition and time.
#' @export
#'
#' @details
#' This function dispatches to either `haz.ExponentialTransition` or `haz.WeibullTransition`
#' based on the `transition` object type.
#'
#' @examples
#' transition <- exponential_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6)
#' haz(transition, 0.4, 2)
haz <- function(transition, t, trans) {
  UseMethod("haz")
}

#' Hazard Function for Exponential Transition Model
#'
#' @param transition (`ExponentialTransition`)\cr
#'   See [exponential_transition()] for details.
#' @param t (`numeric`)\cr time at which hazard is to be computed.
#' @param trans (`integer`)\cr index specifying the transition type.
#'
#' @return Returns the hazard rate (`numeric`) for the exponential transition at the specified time.
#' @export
#'
#' @details
#' Computes the hazard function specifically for an exponential transition model using provided parameters.
#'
#' @examples
#' transition <- exponential_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6)
#' haz.ExponentialTransition(transition, 0.4, 2)
haz.ExponentialTransition <- function(transition, t, trans) {
  # params (in this order): h01, h02, h12
  params <- unlist(transition$hazards)
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
#' haz.WeibullTransition(transition, 0.4, 2)
haz.WeibullTransition <- function(transition, t, trans) {
  # params (in this order): h01, h02, h12, p01, p02, p12
  params <- c(unlist(transition$hazards), unlist(transition$weibull_rates))
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
#' based on the `transition` object type.
#'
#' @examples
#' transition <- exponential_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6)
#' survTrans(transition, 0.4, 2)
survTrans <- function(transition, t, trans) {
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
#' survTrans.ExponentialTransition(transition, 0.4, 2)
survTrans.ExponentialTransition <- function(transition, t, trans) {
  # params (in this order): h01, h02, h12
  params <- unlist(transition$hazards)
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
#' survTrans.WeibullTransition(transition, 0.4, 2)
survTrans.WeibullTransition <- function(transition, t, trans) {
  # params (in this order): h01, h02, h12, p01, p02, p12
  params <- c(unlist(transition$hazards), unlist(transition$weibull_rates))
  exp(-params[trans] * t^params[trans + 3])
}

getInitial <- function(transition) {
  UseMethod("getInitial")
}

getInitial.ExponentialTransition <- function(transition) {
  unlist(transition$hazards)
}

getInitial.WeibullTransition <- function(transition) {
  c(unlist(transition$hazards), unlist(transition$weibull_rates))
}

getTarget <- function(params, data, transition) {
  if ("ExponentialTransition" %in% class(transition)) {
    negLogLik(transition = exponential_transition(h01 = params[1], h02 = params[2], h12 = params[3]), data = data)
  } else {
    negLogLik(transition = weibull_transition(
      h01 = params[1], h02 = params[2], h12 = params[3],
      p01 = params[4], p02 = params[5], p12 = params[6]
    ), data = data)
  }
}

getResults <- function(transition, res) {
  UseMethod("getResults")
}

getResults.ExponentialTransition <- function(transition, res) {
  exponential_transition(h01 = res[1], h02 = res[2], h12 = res[3])
}

getResults.WeibullTransition <- function(transition, res) {
  weibull_transition(
    h01 = res[1], h02 = res[2], h12 = res[3],
    p01 = res[4], p02 = res[5], p12 = res[6]
  )
}

estimateParams <- function(data, transition) {
  data <- prepareData(data)

  res <- stats::optim(
    par = getInitial(transition),
    fn = getTarget,
    method = "L-BFGS-B",
    lower = 1e-3,
    data = data,
    transition = transition
  )$par

  getResults(transition, res)
}
