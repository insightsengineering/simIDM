#' Time-point by which a specified number of events occurred.
#'
#' This returns the study time-point by which a specified number of events (OS or PFS) occurred.
#'
#' @param data  (`data.frame`)\cr illness-death data set in `1rowPatient` format.
#' @param eventNum (`int`)\cr number of events.
#' @param typeEvent (`string`)\cr type of event. Possible values are `OS` and `PFS`.
#' @param byArm  (`logical`)\cr if `TRUE` time-point per treatment arm, else joint evaluation of treatment arms.
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
  if (byArm == FALSE) {
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

  timePoint <- sapply(byVar, function(trt) {
    dataTemp <- data[data$trt == trt, ]
    sortedTimes <- sort(dataTemp$time[dataTemp$event == 1])
    if (eventNum > length(sortedTimes)) {
      cat("Less events occur till EOS than specified \n")
    }
    timePoint <- sortedTimes[eventNum]
  })
  return(timePoint)
}

#' Event-driven censoring.
#'
#' This function censors a study after a pre-specified number of events occurred.
#'
#' @param data (`data.frame`)\cr illness-death data set in `1rowPatient` format.
#' @param eventNum  (`int`)\cr number of events.
#' @param typeEvent (`string`)\cr type of event. Possible values are `OS` and `PFS`.
#'
#' @return This function returns a data set that is censored after `eventNum` of `typeEvent`-events occurred.
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
#' censoringByNumberEvents(data = simStudyWide, eventNum = 20, typeEvent = "PFS")
censoringByNumberEvents <- function(data, eventNum, typeEvent) {
  # first step get timePoint (calendar time) for censoring.
  censTime <- getTimePoint(data, eventNum, typeEvent, byArm = FALSE)
  # censoring time at individual time scale.
  data$censTimeInd <- censTime - data$recruitTime

  # PFStime censored?
  PFScensored <- data$id[data$PFStime > data$censTimeInd]
  # OStime censored?
  OScensored <- data$id[data$OStime > data$censTimeInd]

  # patients that are not yet recruited at censoring time has to be deleted.
  NotRecruited <- data$id[data$censTimeInd < 0]
  censoredData <- data[!(data$id %in% NotRecruited), ]

  # keep minimum of censoring and event time.
  censoredData$OStime <- pmin(censoredData$OStime, censoredData$censTimeInd)
  censoredData$PFStime <- pmin(censoredData$PFStime, censoredData$censTimeInd)

  # adjust event and censoring indicators.
  censoredData$OSevent <- ifelse(censoredData$id %in% OScensored, 0, censoredData$OSevent)
  censoredData$PFSevent <- ifelse(censoredData$id %in% PFScensored, 0, censoredData$PFSevent)

  censoredData$CensoredOS <- abs(1 - censoredData$OSevent)
  censoredData$CensoredPFS <- abs(1 - censoredData$PFSevent)

  # calculate corresponding calendar times.
  censoredData$OStimeCal <- censoredData$OStime + censoredData$recruitTime
  censoredData$PFStimeCal <- censoredData$PFStime + censoredData$recruitTime

  return(censoredData)
}



#' Event tracking in an oncology trial.
#'
#' @param data (`data.frame`)\cr illness-death data set in `1rowPatient` format.
#' @param timeP (`numeric`)\cr  study time-point.
#' @param byArm  (`logical`)\cr  if `TRUE` time-point per treatment arm, else joint evaluation of treatment arms.
#'
#' @return This function returns a data frame including number of PFS events, number of OS events,
#'  number of recruited patients, number of censored patients and number of ongoing patients at `timeP`.
#'
#'
#' @export
#'
#' @examples
#' simStudy <- getOneClinicalTrial(
#'   nPat = c(20, 20), transitionByArm = list(transition1, transition2),
#'   dropout = list(rate = 0.3, time = 10),
#'   accrual = list(param = "time", value = 0)
#' )
#' simStudyWide <- getDatasetWideFormat(simStudy)
#' trackEventsPerTrial(data = simStudyWide, timeP = 1.5, byArm = FALSE)
trackEventsPerTrial <- function(data, timeP, byArm = FALSE) {
  if (byArm == FALSE) {
    data$trt <- "all"
  }
  byVar <- unique(data$trt)
  allNumbers <- lapply(byVar, function(j) {
    datTemp <- data[data$trt == j, ]
    eventsPFS <- sapply(timeP, function(t) {
      return(sum(datTemp$PFSevent[(datTemp$PFStime + datTemp$recruitTime) <= t]))
    })
    eventsOS <- sapply(timeP, function(t) {
      return(sum(datTemp$OSevent[(datTemp$OStime + datTemp$recruitTime) <= t]))
    })
    eventsTrial <- sapply(timeP, function(t) {
      allRecruited <- datTemp$id[datTemp$recruitTime <= t]
      allCensored <- datTemp$id[(datTemp$OStime + datTemp$recruitTime) <= t & datTemp$OSevent == 0]
      allDeath <- datTemp$id[(datTemp$OStime + datTemp$recruitTime) <= t & datTemp$OSevent == 1]
      ## recruited, not death, not censored.
      allUnderObs <- datTemp$id[datTemp$id %in% allRecruited & !(datTemp$id %in% allCensored) &
        !(datTemp$id %in% allDeath)]
      numRecruited <- length(unique(allRecruited))
      numCensored <- length(unique(allCensored))
      numUnderObs <- length(unique(allUnderObs))

      return(c(Recruited = numRecruited, Censored = numCensored, UnderObs = numUnderObs))
    })


    allNumbers <- rbind(eventsPFS, eventsOS, eventsTrial)
    rownames(allNumbers) <- c("PFS", "OS", "Recruited", "Censored", "Ongoing")
    colnames(allNumbers) <- paste0("Timepoint: ", round(timeP, 2))
    return(allNumbers)
  })
  names(allNumbers) <- byVar

  return(allNumbers)
}
