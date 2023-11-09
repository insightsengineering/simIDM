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


# transition <- weibull_transition(h01 = 1.5, h02 = 0.5, h12 = 1.1, p01 = 0.4, p02 = 0.5, p12 = 2.4)
# trial <- getOneClinicalTrial(nPat = c(10000),
#                             transitionByArm = list(transition),
#                             dropout = list(rate = 0.3, time = 1),
#                             accrual = list(param = "intensity", value = 100))

# estimateParams(trial, "Weibull")

negLogLik <- function(transition, data) {
  with(data, -sum(log(haz(transition, exit, trans)^status * survTrans(transition, exit, trans) / survTrans(transition, entry, trans))))
}

haz <- function(transition, t, trans = c(1, 2, 3)) {
  trans <- match.arg(trans)
  UseMethod(haz)
}

haz.ExponentialTransition <- function(transition, t, trans = c(1, 2, 3)) {
  # params (in this order): h01, h02, h12
  params <- unlist(transition$hazards)
  exp(-params[trans] * t)
}

haz.WeibullTransition <- function(transition, t, trans = c(1, 2, 3)) {
  # params (in this order): h01, h02, h12, p01, p02, p12
  params <- c(unlist(transition$hazards), unlist(transition$weibull_rates))
  exp(-params[trans] * t)
}

survTrans <- function(transition, t, trans = c(1, 2, 3)) {
  trans <- match.arg(trans)
  UseMethod(survTrans)
}

survTrans.Exponential <- function(transition, t, trans = c(1, 2, 3)) {
  # params (in this order): h01, h02, h12
  params <- unlist(transition$hazards)
  exp(-params[trans] * t)
}

survTrans.Weibull <- function(transition, t, trans = c(1, 2, 3)) {
  # params (in this order): h01, h02, h12, p01, p02, p12
  params <- c(unlist(transition$hazards), unlist(transition$weibull_rates))
  exp(-params[trans] * t^params[trans + 3])
}

getNegLogLik <- function(params, data, family = c("exponential", "Weibull")) {
  family <- match.arg(family)
  if (family == "exponential") {
    negLogLik(transition = exponential_transition(h01 = params[1], h02 = params[2], h12 = params[3]), data = data)
  } else {
    negLogLik(transition = weibull_transition(
      h01 = params[1], h02 = params[2], h12 = params[3],
      p01 = params[4], p02 = params[5], p12 = params[6]
    ), data = data)
  }
}

estimateParams <- function(data, family = c("exponential", "Weibull")) {
  data <- prepareData(data)
  family <- match.arg(family)
  initial <- if (family == "exponential") {
    c(1, 1, 1)
  } else {
    c(1, 1, 1, 1, 1, 1)
  }

  res <- stats::optim(
    par = initial,
    getNegLogLik,
    method = "Nelder-Mead",
    data = data,
    family = family
  )$par

  if (family == "exponential") {
    list("h01" = res[1], "h02" = res[2], "h12" = res[3])
  } else {
    list(
      "h01" = res[1], "h02" = res[2], "h12" = res[3],
      "p01" = res[4], "p02" = res[5], "p12" = res[6]
    )
  }
}
