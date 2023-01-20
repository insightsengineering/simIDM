#' Simulation of a Single Oncology Clinical Trial
#'
#' This function creates a data set with a single simulated oncology clinical trial with one row per transition
#' based on an illness-death model. Studies with an arbitrary number of treatment arms are possible.
#'
#' @param nPat (`integer`)\cr numbers of patients per treatment arm.
#' @param transitionByArm (`list`) \cr transition parameters for each treatment group.
#'   See [exponential_transition()], [piecewise_exponential()] and [weibull_transition()] for details.
#' @param dropout  dropout (`list`)\cr specifies drop-out probability. See [getSimulatedData()] for details.
#'  Can be specified either as one list that should be applied to all treatment groups or a separate list
#'  for each treatment group.
#' @param accrual  accrual (`list`)\cr specifies accrual intensity. See [addStaggeredEntry()] for details.
#'  Can be specified either as one list that should be applied to all treatment groups or a separate list
#'  for each treatment group.
#'
#' @return This returns a data frame with one simulated clinical trial and multiple treatment arms.
#'   See [getSimulatedData()] for the explanation of the columns. The column `trt` contains the treatment indicator.
#'   This is a helper function of [getClinicalTrials()].
#' @export
#'
#' @examples
#' transition1 <- exponential_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6)
#' transition2 <- exponential_transition(h01 = 1, h02 = 1.3, h12 = 1.7)
#' transition3 <- exponential_transition(h01 = 1.1, h02 = 1, h12 = 1.5)
#' getOneClinicalTrial(
#'   nPat = c(30, 20, 30), transitionByArm = list(transition1, transition2, transition3),
#'   dropout = list(rate = 0, time = 12),
#'   accrual = list(param = "time", value = 0)
#' )
getOneClinicalTrial <- function(nPat, transitionByArm,
                                dropout = list(rate = 0, time = 12),
                                accrual = list(param = "time", value = 0)) {
  assert_list(transitionByArm)
  nPat <- as.integer(nPat)
  nArm <- length(transitionByArm)
  assert_integer(nPat,
    lower = 1,
    any.missing = FALSE,
    all.missing = FALSE,
    len = nArm,
  )
  assert_list(dropout)
  assert_list(accrual)

  # Same accrual and dropout parameters for each group?
  if (is.list(dropout[[1]])) {
    assert_list(dropout, len = nArm, types = "list")
  } else {
    dropout <- rep(list(dropout), nArm)
  }

  if (is.list(accrual[[1]])) {
    assert_list(accrual, len = nArm, types = "list")
  } else {
    accrual <- rep(list(accrual), nArm)
  }

  # Each loop simulates a single trial arm.
  # Starting values for the loop.
  simdata <- NULL
  previousPts <- 0
  for (i in seq_len(nArm)) {
    group <- getSimulatedData(nPat[i], transitionByArm[[i]], dropout[[i]], accrual[[i]])
    group$trt <- i
    group$id <- group$id + previousPts
    simdata <- rbind(simdata, group)
    previousPts <- previousPts + nPat[i]
  }
  simdata
}

#' Conversion of a Data Set from One Row per Transition to One Row per Patient
#'
#' @param data (`data.frame`)\cr data frame containing entry and exit times of an illness-death model.
#'   See [getSimulatedData()] for details.
#'
#' @return This function returns a data set with one row per patient and endpoints PFS and OS.
#' @export
#'
#' @details
#' The output data set contains the following columns:
#' - id (`integer`): patient id.
#' - trt `integer`): treatment id.
#' - PFStime (`numeric`): event time of PFS event.
#' - CensoredPFS (`logical`): censoring indicator for PFS event.
#' - PFSevent (`logical`): event indicator for PFS event.
#' - OStime (`numeric`): event time of OS event.
#' - CensoredOS (`logical`): censoring indicator for OS event.
#' - OSevent (`logical`): event indicator for OS event.
#' - recruitTime (`numeric`): time of recruitment.
#' - OStimeCal (`numeric`): OS event time at calendar time scale.
#' - PFStimeCal (`numeric`): PFS event time at calendar time scale.
#'
#' @examples
#' transition1 <- exponential_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6)
#' transition2 <- exponential_transition(h01 = 1, h02 = 1.3, h12 = 1.7)
#' transition3 <- exponential_transition(h01 = 1.1, h02 = 1, h12 = 1.5)
#' simData <- getOneClinicalTrial(
#'   nPat = c(30, 20, 30), transitionByArm = list(transition1, transition2, transition3),
#'   dropout = list(rate = 0, time = 12),
#'   accrual = list(param = "time", value = 0)
#' )
#' getDatasetWideFormat(simData)
getDatasetWideFormat <- function(data) {
  assert_data_frame(data, ncols = 9)
  assert_subset(c("id", "from", "to", "entry", "exit", "entryAct", "exitAct", "censAct", "trt"), names(data))

  # Recruitment time is the actual entry time of the initial state.
  recruitTime <- subset(data[, c("id", "entryAct")], data$from == 0)
  names(recruitTime)[names(recruitTime) == "entryAct"] <- "recruitTime"

  # The OS time is the entry time into state 2 or the censoring time whatever occurs first.
  OStime <- subset(data[, c("id", "exit")], data$to == 2 | data$to == "cens")
  names(OStime)[names(OStime) == "exit"] <- "OStime"

  # The PFS time is the entry time into state 1 or state 2 or the censoring time whatever occurs first.
  PFStime <- subset(
    data[, c("id", "exit")],
    (data$to == 2 & data$from == 0) | (data$to == 1) | (data$to == "cens" & data$from == 0)
  )
  names(PFStime)[names(PFStime) == "exit"] <- "PFStime"

  # Add censoring indicators.
  censoredIdsPFS <- data$id[data$to == "cens" & data$from == 0]
  censoredIdsOS <- data$id[data$to == "cens"]
  id <- unique(data$id)
  CensoredOS <- cbind(id = id, CensoredOS = as.integer(id %in% censoredIdsOS))
  CensoredPFS <- cbind(id = id, CensoredPFS = as.integer(id %in% censoredIdsPFS))

  # Merge all data sets to one.
  newdata <- unique(data[, c("id", "trt")])
  newdata <- merge(x = newdata, y = PFStime, by = "id")
  newdata <- merge(x = newdata, y = CensoredPFS, by = "id")

  # Do we have an observed PFS event?
  newdata$PFSevent <- abs(1 - newdata$CensoredPFS)
  newdata <- merge(x = newdata, y = OStime, by = "id")
  newdata <- merge(x = newdata, y = CensoredOS, by = "id")

  # Do we have an observed OS event?
  newdata$OSevent <- abs(1 - newdata$CensoredOS)
  newdata <- merge(x = newdata, y = recruitTime, by = "id")

  # Add variables with event times in calendar time.
  newdata$OStimeCal <- newdata$OStime + newdata$recruitTime
  newdata$PFStimeCal <- newdata$PFStime + newdata$recruitTime

  newdata
}

#' Simulation of a Large Number of Oncology Clinical Trials
#'
#' @param nRep (`int`)\cr number of simulated trials.
#' @param ... parameters transferred to [getOneClinicalTrial()], see  [getOneClinicalTrial()] for details.
#' @param seed (`int`)\cr random seed used for this simulation.
#' @param datType (`string`)\cr possible values are `1rowTransition` and `1rowPatient`.
#'
#' @return This function returns a list with `nRep` simulated data sets in the format specified by `datType`.
#'  See [getDatasetWideFormat()] [getOneClinicalTrial()] for details.
#' @export
#'
#' @examples
#' transition1 <- exponential_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6)
#' transition2 <- exponential_transition(h01 = 1, h02 = 1.3, h12 = 1.7)
#' getClinicalTrials(
#'   nRep = 10, nPat = c(20, 20), seed = 1234, datType = "1rowTransition",
#'   transitionByArm = list(transition1, transition2), dropout = list(rate = 0.5, time = 12),
#'   accrual = list(param = "intensity", value = 7)
#' )
getClinicalTrials <- function(nRep, ..., seed = 1234, datType = "1rowTransition") {
  assert_number(nRep, lower = 1)
  assert_choice(datType, c("1rowTransition", "1rowPatient"))

  set.seed(seed)
  # getOneClinicalTrial generates a single clinical trial with multiple arms. Generate nRep simulated trials:
  simulatedTrials <- lapply(
    seq_len(nRep),
    FUN = function(x, ...) getOneClinicalTrial(...),
    ...
  )

  # Final data set format: one row per patient or one row per transition?
  if (datType == "1rowPatient") {
    simulatedTrials <- lapply(simulatedTrials, getDatasetWideFormat)
  }
  simulatedTrials
}
