#' Time-point by which a specified number of events occurred.
#'
#' This returns the study time-point by which a specified number of events (PFS or OS) occurred.
#'
#' @param data (`data.frame`)\cr illness-death data set in `1rowPatient` format.
#' @param eventNum (`int`)\cr number of events.
#' @param typeEvent (`string`)\cr type of event. Possible values are `PFS` and `OS`.
#' @param byArm  (`logical`)\cr if `TRUE` time-point per treatment arm, else joint evaluation
#'   of treatment arms.
#'
#' @return This returns  the time-point by which `eventNum` of `typeEvent`-events occurred.
#' @export
#'
#' @examples
#' transition1 <- weibull_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6, p01 = 0.8, p02 = 0.9, p12 = 1)
#' transition2 <- weibull_transition(h01 = 1, h02 = 1.3, h12 = 1.7, p01 = 1.1, p02 = 0.9, p12 = 1.1)
#'
#' simStudy <- getOneClinicalTrial(
#'   nPat = c(20, 20), transitionByArm = list(transition1, transition2),
#'   dropout = list(rate = 0.3, time = 10),
#'   accrual = list(param = "time", value = 0)
#' )
#' simStudyWide <- getDatasetWideFormat(simStudy)
#' getTimePoint(simStudyWide, eventNum = 10, typeEvent = "OS", byArm = FALSE)
getTimePoint <- function(data, eventNum, typeEvent, byArm = FALSE) {
  assert_data_frame(data, ncols = 11)
  assert_int(eventNum, lower = 1L)
  assert_choice(typeEvent, c("OS", "PFS"))
  assert_logical(byArm)

  if (!byArm) {
    data$trt <- "all"
  }
  byVar <- unique(data$trt)

  if (typeEvent == "OS") {
    data$event <- data$OSevent
    data$time <- data$OStimeCal
  } else if (typeEvent == "PFS") {
    data$event <- data$PFSevent
    data$time <- data$PFStimeCal
  }

  sapply(byVar, function(trt) {
    dataTemp <- data[data$trt == trt, ]
    sortedTimes <- sort(dataTemp$time[dataTemp$event == 1])
    if (eventNum > length(sortedTimes)) {
      warning("Less events occur till EOS than specified \n")
    }
    timePoint <- sortedTimes[eventNum]
  })
}

#' Helper function for `censoringByNumberEvents`
#'
#' @param time (`numeric`) \cr event times.
#' @param event (`numeric`)\cr event indicator.
#' @param data (`data.frame`)\cr data frame including patient id `id`, recruiting time `recruitTime`
#'  and individual censoring time `censTimeInd`.
#'
#' @return This function returns a data frame with columns:
#' event time, censoring indicator, event indicator and event time
#' in calendar time.
#' @export
#'
#' @examples
#' transition1 <- weibull_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6, p01 = 0.8, p02 = 0.9, p12 = 1)
#' transition2 <- weibull_transition(h01 = 1, h02 = 1.3, h12 = 1.7, p01 = 1.1, p02 = 0.9, p12 = 1.1)
#'
#' simStudy <- getOneClinicalTrial(
#'   nPat = c(20, 20), transitionByArm = list(transition1, transition2),
#'   dropout = list(rate = 0.3, time = 10),
#'   accrual = list(param = "time", value = 7)
#' )
#' simStudyWide <- getDatasetWideFormat(simStudy)
#' simStudyWide$censTimeInd <- 5 - simStudyWide$recruitTime
#' NotRecruited <- simStudyWide$id[simStudyWide$censTimeInd < 0]
#' censoredData <- simStudyWide[!(simStudyWide$id %in% NotRecruited), ]
#' getCensoredData(time = censoredData$OStime, event = censoredData$OSevent, data = censoredData)
getCensoredData <- function(time, event, data) {
  assert_numeric(time, lower = 0)
  assert_numeric(event, lower = 0, upper = 1)
  assert_data_frame(data)

  # Event time censored?
  epsilon <- 1e-10
  Censored <- data$id[(time - epsilon) > data$censTimeInd]
  # keep minimum of censoring and event time. epsilon)
  Censoredtime <- pmin(time, data$censTimeInd)
  # adjust event and censoring indicators.
  CensoredEvent <- ifelse(data$id %in% Censored, 0, event)
  Censored <- abs(1 - CensoredEvent)
  # calculate corresponding calendar times.
  timeCal <- Censoredtime + data$recruitTime

  data.frame(
    time = Censoredtime,
    Censored = Censored,
    event = CensoredEvent,
    timeCal = timeCal
  )
}

#' Event-driven censoring.
#'
#' This function censors a study after a pre-specified number of events occurred.
#'
#' @param data (`data.frame`)\cr illness-death data set in `1rowPatient` format.
#' @param eventNum  (`int`)\cr number of events.
#' @param typeEvent (`string`)\cr type of event. Possible values are `PFS` and `OS`.
#'
#' @return This function returns a data set that is censored after `eventNum` of
#'   `typeEvent`-events occurred.
#' @export
#'
#' @examples
#' transition1 <- weibull_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6, p01 = 0.8, p02 = 0.9, p12 = 1)
#' transition2 <- weibull_transition(h01 = 1, h02 = 1.3, h12 = 1.7, p01 = 1.1, p02 = 0.9, p12 = 1.1)
#'
#' simStudy <- getOneClinicalTrial(
#'   nPat = c(20, 20), transitionByArm = list(transition1, transition2),
#'   dropout = list(rate = 0.3, time = 10),
#'   accrual = list(param = "time", value = 7)
#' )
#' simStudyWide <- getDatasetWideFormat(simStudy)
#' censoringByNumberEvents(data = simStudyWide, eventNum = 20, typeEvent = "PFS")
censoringByNumberEvents <- function(data, eventNum, typeEvent) {
  assert_data_frame(data, ncols = 11)
  assert_int(eventNum, lower = 1L)
  assert_choice(typeEvent, c("OS", "PFS"))

  # first step get timePoint (calendar time) for censoring.
  censTime <- getTimePoint(data, eventNum, typeEvent, byArm = FALSE)

  # censoring time at individual time scale.
  data$censTimeInd <- censTime - data$recruitTime
  # patients that are not yet recruited at censoring time has to be deleted.
  NotRecruited <- data$id[data$censTimeInd < 0]
  censoredData <- data[!(data$id %in% NotRecruited), ]

  OSdata <- getCensoredData(time = censoredData$OStime, event = censoredData$OSevent, data = censoredData)
  PFSdata <- getCensoredData(time = censoredData$PFStime, event = censoredData$PFSevent, data = censoredData)

  data.frame(
    id = censoredData$id, trt = censoredData$trt,
    PFStime = PFSdata$time, PFSevent = PFSdata$event,
    OStime = OSdata$time, CensoredOS = OSdata$Censored, OSevent = OSdata$event,
    recruitTime = censoredData$recruitTime, OStimeCal = OSdata$timeCal, PFStimeCal = PFSdata$timeCal
  )
}

#' Number of recruited/censored/ongoing Patients.
#'
#'
#' @param data (`data.frame`)\cr illness-death data set in `1rowPatient` format.
#' @param t (`numeric`)\cr  study time-point.
#'
#' @return This function returns number of recruited patients,
#' number of censored and number of patients under observations.
#' @export
#'
#' @examples
#' transition1 <- weibull_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6, p01 = 0.8, p02 = 0.9, p12 = 1)
#' transition2 <- weibull_transition(h01 = 1, h02 = 1.3, h12 = 1.7, p01 = 1.1, p02 = 0.9, p12 = 1.1)
#'
#' simStudy <- getOneClinicalTrial(
#'   nPat = c(20, 20), transitionByArm = list(transition1, transition2),
#'   dropout = list(rate = 0.6, time = 10),
#'   accrual = list(param = "time", value = 0)
#' )
#' simStudyWide <- getDatasetWideFormat(simStudy)
#' getEventsAll(data = simStudyWide, t = 1.5)
getEventsAll <- function(data, t) {
  assert_data_frame(data, ncols = 11)
  assert_positive_number(t, zero_ok = FALSE)

  allRecruited <- data$id[data$recruitTime <= t]
  allCensored <- data$id[(data$OStime + data$recruitTime) <= t & data$OSevent == 0]
  allDeath <- data$id[(data$OStime + data$recruitTime) <= t & data$OSevent == 1]
  ## recruited, not death, not censored.
  allUnderObs <- data$id[data$id %in% allRecruited & !(data$id %in% allCensored) &
    !(data$id %in% allDeath)]
  numRecruited <- length(unique(allRecruited))
  numCensored <- length(unique(allCensored))
  numUnderObs <- length(unique(allUnderObs))

  c(
    Recruited = numRecruited,
    Censored = numCensored,
    UnderObs = numUnderObs
  )
}

#' Helper Function for `trackEventsPerTrial`
#'
#' @param event (`numeric`)\cr event indicator.
#' @param time (`numeric`) \cr event times.
#' @param t  (`numeric`)\cr  study time-point.
#'
#' @return This function returns the number of events occurred until time t.
#' @export
#'
#' @examples
#' event <- c(0, 1, 1, 1, 0)
#' time <- c(3, 3.4, 5, 6, 5.5)
#' getNumberEvents(event = event, time = time, t = 5)
getNumberEvents <- function(event, time, t) {
  assert_numeric(event, lower = 0, upper = 1)
  assert_numeric(time, lower = 0)
  assert_number(t, lower = 0)

  sum(event[time <= t])
}


#' Event tracking in an oncology trial.
#'
#' @param data (`data.frame`)\cr illness-death data set in `1rowPatient` format.
#' @param timeP (`numeric`)\cr  vector of study time-points.
#' @param byArm  (`logical`)\cr  if `TRUE` time-point per treatment arm, else joint evaluation of treatment arms.
#'
#' @return This function returns a data frame including number of PFS events, number of OS events,
#'  number of recruited patients, number of censored patients and number of ongoing patients at `timeP`.
#'
#' @export
#'
#' @examples
#' transition1 <- weibull_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6, p01 = 0.8, p02 = 0.9, p12 = 1)
#' transition2 <- weibull_transition(h01 = 1, h02 = 1.3, h12 = 1.7, p01 = 1.1, p02 = 0.9, p12 = 1.1)
#'
#' simStudy <- getOneClinicalTrial(
#'   nPat = c(20, 20), transitionByArm = list(transition1, transition2),
#'   dropout = list(rate = 0.3, time = 10),
#'   accrual = list(param = "time", value = 0)
#' )
#' simStudyWide <- getDatasetWideFormat(simStudy)
#' trackEventsPerTrial(data = simStudyWide, timeP = 1.5, byArm = FALSE)
trackEventsPerTrial <- function(data, timeP, byArm = FALSE) {
  assert_data_frame(data, ncols = 11)
  assert_numeric(timeP, lower = 0)
  assert_logical(byArm)

  if (!byArm) {
    data$trt <- "all"
  }
  byVar <- unique(data$trt)
  allNumbers <- lapply(byVar, function(j) {
    datTemp <- data[data$trt == j, ]
    eventsPFS <- sapply(timeP, getNumberEvents, event = datTemp$PFSevent, time = datTemp$PFStimeCal)
    eventsOS <- sapply(timeP, getNumberEvents, event = datTemp$OSevent, time = datTemp$OStimeCal)
    eventsTrial <- sapply(timeP, getEventsAll, data = datTemp)

    allNumbers <- rbind(eventsPFS, eventsOS, eventsTrial)
    rownames(allNumbers) <- c("PFS", "OS", "Recruited", "Censored", "Ongoing")
    colnames(allNumbers) <- paste0("Timepoint: ", round(timeP, 2))
    allNumbers
  })
  names(allNumbers) <- byVar
  allNumbers
}
