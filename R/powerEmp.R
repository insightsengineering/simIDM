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
logRankTest <- function(data, typeEvent, critical) {
  # must be `1rowTransition`, i.e. have 10 (censored) or 11 columns
  assert_data_frame(data, min.cols = 10, max.cols = 11)
  assert_choice(typeEvent, c("OS", "PFS"))
  assert_positive_number(critical)

  if (typeEvent == "OS") {
    time <- data$OStime
    event <- data$OSevent
  } else if (typeEvent == "PFS") {
    time <- data$PFStime
    event <- data$PFSevent
  }

  logRank <- survdiff(Surv(time, event) ~ trt, data)
  sqrt(logRank$chisq) > critical
}

#' Empirical Power for a List of Simulated Trials
#'
#' This function computes four types of empirical power — PFS, OS, at-least (significant
#' in at least one of PFS/OS), and joint (significant in both PFS and OS) — using
#' the log-rank test. Empirical power is calculated as the proportion of significant
#' results in simulated trials, each ending when a set number of PFS/OS events occur.
#' Critical values for PFS and OS test significance must be specified.
#'
#' @param simulatedTrials (`list`)\cr simulated trial data sets, see [getClinicalTrials()].
#' @param criticalValuePFS (positive `number`)\cr critical value of the log-rank test for PFS.
#' @param criticalValueOS (positive `number`)\cr critical value of the log-rank test for OS.
#' @param eventNumPFS (`integer`)\cr number of PFS events required to trigger PFS analysis.
#' @param eventNumOS (`integer`)\cr number of OS events required to trigger OS analysis.
#'
#' @return This returns values of four measures of empirical power.
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
#' powerEmp(
#'   simTrials = simTrials, criticalPFS = 2.4, criticalOS = 2.2,
#'   eventNumPFS = 300, eventNumOS = 500
#' )
powerEmp <- function(simTrials, criticalPFS, criticalOS, eventNumPFS, eventNumOS) {
  assert_list(simTrials, null.ok = FALSE)
  assert_positive_number(criticalPFS)
  assert_positive_number(criticalOS)
  assert_int(eventNumPFS, lower = 1L)
  assert_int(eventNumOS, lower = 1L)

  nRep <- length(simTrials)
  simTrials <- lapply(simTrials, function(x) if (ncol(x) == 9) getDatasetWideFormat(x) else x)

  # Censor simulated trials at time-point of OS/PFS analysis
  trialsAnaPFS <- lapply(simTrials, censoringByNumberEvents,
    eventNum = eventNumPFS, typeEvent = "PFS"
  )
  trialsAnaOS <- lapply(simTrials, censoringByNumberEvents,
    eventNum = eventNumOS, typeEvent = "OS"
  )

  # Compute log-rank test for all trials for both OS/PFS
  passedLogRankPFS <- unlist(lapply(trialsAnaPFS,
    logRankTest,
    typeEvent = "PFS",
    criticalPFS
  ))
  passedLogRankOS <- unlist(lapply(trialsAnaOS,
    logRankTest,
    typeEvent = "OS",
    criticalOS
  ))

  # Empirical power is the fraction of trials with significant log-rank test
  powerPFS <- sum(passedLogRankPFS) / nRep
  powerOS <- sum(passedLogRankOS) / nRep

  # Derived measures of power
  sumPassed <- passedLogRankPFS + passedLogRankOS
  # At-least power: At least one of PFS/OS has significant log-rank test
  powerAtLeast <- sum(sumPassed >= 1) / nRep
  # Joint power: Both PFS/OS have significant log-rank tests
  powerJoint <- sum(sumPassed == 2) / nRep

  c("powerPFS" = powerPFS, "powerOS" = powerOS, "powerAtLeast" = powerAtLeast, "powerJoint" = powerJoint)
}
