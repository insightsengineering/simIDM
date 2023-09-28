#' Log-Rank Test for a Single Trial
#'
#' This function evaluates the significance of either PFS or OS endpoints in a trial,
#' based on a pre-specified critical value.
#'
#' @param data (`data.frame`)\cr data frame containing entry and exit times of an
#' illness-death model. See [getSimulatedData()] for details.
#' @param typeEvent (`string`)\cr endpoint to be evaluated, possible values are `PFS` and `OS`.
#' @param critical (positive `number`)\cr critical value of the log-rank test.
#'
#' @return Logical value indicating log-rank test significance.
#' @export
#'
#' @examples
#' transition1 <- exponential_transition(h01 = 0.06, h02 = 0.3, h12 = 0.3)
#' transition2 <- exponential_transition(h01 = 0.1, h02 = 0.4, h12 = 0.3)
#' simTrial <- getClinicalTrials(
#'   nRep = 1, nPat = c(800, 800), seed = 1234, datType = "1rowPatient",
#'   transitionByArm = list(transition1, transition2), dropout = list(rate = 0.5, time = 12),
#'   accrual = list(param = "intensity", value = 7)
#' )[[1]]
#' logRankTest(data = simTrial, typeEvent = "OS", critical = 3.4)
logRankTest <- function(data, typeEvent = c("PFS", "OS"), critical) {
  # Must be have the format of one row per patient (`datType` must be 1rowPatient`
  # in getClinicalTrials()), i.e. have 10 (if censored) or 11 columns.
  assert_data_frame(data, min.cols = 10, max.cols = 11)
  typeEvent <- match.arg(typeEvent)
  assert_positive_number(critical)

  if (typeEvent == "OS") {
    time <- data$OStime
    event <- data$OSevent
  } else if (typeEvent == "PFS") {
    time <- data$PFStime
    event <- data$PFSevent
  }

  logRank <- survival::survdiff(survival::Surv(time, event) ~ trt, data)
  sqrt(logRank$chisq) > critical
}


#' Helper function to conduct log-rank tests for either PFS or OS
#'
#' This function evaluates the significance of either PFS or OS endpoints for each trial
#' in a list of trials, based on a pre-specified critical value.
#'
#' @param simTrials (`list`)\cr simulated trial data sets, see [getClinicalTrials()].
#' @param typeEvent (`string`)\cr endpoint to be evaluated, possible values are `PFS` and `OS`.
#' @param eventNum (`integer`)\cr number of events required to trigger analysis.
#' @param critical (positive `number`)\cr critical value of the log-rank test.
#'
#' @return Logical vector indicating log-rank test significance for each trial.
#' @export
#'
#' @examples
#' transition1 <- exponential_transition(h01 = 0.06, h02 = 0.3, h12 = 0.3)
#' transition2 <- exponential_transition(h01 = 0.1, h02 = 0.4, h12 = 0.3)
#' simTrials <- getClinicalTrials(
#'   nRep = 50, nPat = c(800, 800), seed = 1234, datType = "1rowPatient",
#'   transitionByArm = list(transition1, transition2), dropout = list(rate = 0.5, time = 12),
#'   accrual = list(param = "intensity", value = 7)
#' )
#' passedLogRank(simTrials = simTrials, typeEvent = "PFS", eventNum = 300, critical = 2.4)
#' @keywords internal
passedLogRank <- function(simTrials, typeEvent, eventNum, critical) {
  assert_list(simTrials, null.ok = FALSE)
  assert_positive_number(critical)
  assert_count(eventNum, positive = TRUE)

  # Censor simulated trials at time-point of OS/PFS analysis.
  trialsAna <- lapply(
    X = simTrials,
    FUN = censoringByNumberEvents,
    eventNum = eventNum,
    typeEvent = typeEvent
  )
  # Compute log-rank test for all trials for OS/PFS.
  unlist(lapply(
    X = trialsAna,
    FUN = logRankTest,
    typeEvent = typeEvent,
    critical = critical
  ))
}

#' Empirical Significance for a List of Simulated Trials
#'
#' This function computes four types of empirical significance — PFS, OS, at-least (significant
#' in at least one of PFS/OS), and joint (significant in both PFS and OS) — using
#' the log-rank test. Empirical significance is calculated as the proportion of significant
#' results in simulated trials, each ending when a set number of PFS/OS events occur.
#' Critical values for PFS and OS test significance must be specified. If trials simulate equal
#' transition hazards across groups (H0), empirical significance estimates type I error;
#' if they simulate differing transition hazards (H1), it estimates power.
#'
#' @param simTrials (`list`)\cr simulated trial data sets, see [getClinicalTrials()].
#' @param criticalPFS (positive `number`)\cr critical value of the log-rank test for PFS.
#' @param criticalOS (positive `number`)\cr critical value of the log-rank test for OS.
#' @param eventNumPFS (`integer`)\cr number of PFS events required to trigger PFS analysis.
#' @param eventNumOS (`integer`)\cr number of OS events required to trigger OS analysis.
#'
#' @return This returns values of four measures of empirical significance.
#' @export
#'
#' @examples
#' transition1 <- exponential_transition(h01 = 0.06, h02 = 0.3, h12 = 0.3)
#' transition2 <- exponential_transition(h01 = 0.1, h02 = 0.4, h12 = 0.3)
#' simTrials <- getClinicalTrials(
#'   nRep = 50, nPat = c(800, 800), seed = 1234, datType = "1rowPatient",
#'   transitionByArm = list(transition1, transition2), dropout = list(rate = 0.5, time = 12),
#'   accrual = list(param = "intensity", value = 7)
#' )
#' empSignificant(
#'   simTrials = simTrials, criticalPFS = 2.4, criticalOS = 2.2,
#'   eventNumPFS = 300, eventNumOS = 500
#' )
empSignificant <- function(simTrials, criticalPFS, criticalOS, eventNumPFS, eventNumOS) {
  assert_list(simTrials, null.ok = FALSE)
  assert_positive_number(criticalPFS)
  assert_positive_number(criticalOS)
  assert_count(eventNumPFS, positive = TRUE)
  assert_count(eventNumOS, positive = TRUE)

  nRep <- length(simTrials)
  simTrials <- lapply(
    X = simTrials,
    FUN = function(x) if (ncol(x) == 9) getDatasetWideFormat(x) else x
  )

  # Which trials passed the log-rank test for PFS/OS?
  passedLogRankPFS <- passedLogRank(
    simTrials = simTrials,
    typeEvent = "PFS",
    eventNum = eventNumPFS,
    critical = criticalPFS
  )
  passedLogRankOS <- passedLogRank(
    simTrials = simTrials,
    typeEvent = "OS",
    eventNum = eventNumOS,
    critical = criticalOS
  )

  # Empirical significance is the fraction of trials with significant log-rank test:
  significantPFS <- sum(passedLogRankPFS) / nRep
  significantOS <- sum(passedLogRankOS) / nRep

  # Derived measures of significance:
  sumPassed <- passedLogRankPFS + passedLogRankOS
  # At least one of PFS/OS has significant log-rank test.
  significantAtLeastOne <- sum(sumPassed >= 1) / nRep
  # Both PFS/OS have significant log-rank tests.
  significantBoth <- sum(sumPassed == 2) / nRep

  list(
    "significantPFS" = significantPFS,
    "significantOS" = significantOS,
    "significantAtLeastOne" = significantAtLeastOne,
    "significantBoth" = significantBoth
  )
}
