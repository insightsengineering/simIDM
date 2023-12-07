#' Transitions from the Intermediate State to the Absorbing State
#'
#' This function creates transition entry and exit times from the intermediate state to the absorbing state
#' for an existing data frame containing the exit times out of the initial state.
#'
#' @param simDataOne (`data.frame`)\cr a data frame containing all patients with transitions
#'  into the intermediate state. See [getSimulatedData()] for details.
#' @param transition (`TransitionParameters`)\cr transition parameters comprising
#'  `hazards`, corresponding `intervals` and `weibull_rates`, see [exponential_transition()], [piecewise_exponential()]
#'  and [weibull_transition()] for details.
#'
#' @return This returns a data frame with one row per patient for the second transition,
#'  i.e. the transition out of the intermediate
#'  state. This is a helper function of [getSimulatedData()].
#'
#' @export
#'
#' @examples
#' simDataOne <- data.frame(
#'   id = c(1:3), to = c(1, 1, 1), from = c(0, 0, 0), entry = c(0, 0, 0),
#'   exit = c(3, 5.6, 7.2), censTime = c(6.8, 5.9, 9.4)
#' )
#' transition <- exponential_transition(1, 1.6, 0.3)
#' getOneToTwoRows(simDataOne, transition)
getOneToTwoRows <- function(simDataOne, transition) {
  assert_data_frame(simDataOne, ncols = 6)
  assert_class(transition, "TransitionParameters")

  id1 <- simDataOne$id
  N1 <- nrow(simDataOne)
  U1 <- stats::runif(N1)
  to1 <- rep(2, N1)
  from1 <- rep(1, N1)
  entry1 <- simDataOne$exit
  h12 <- transition$hazard$h12
  p12 <- transition$weibull_rates$p12
  pw12 <- transition$intervals$pw12

  # Create waiting time in state 1.
  wait_time1 <- if (transition$family == "exponential") {
    (-log(1 - U1)) / (h12)
  } else if (transition$family == "Weibull") {
    getWaitTimeSum(U = U1, haz1 = h12, haz2 = 0, p1 = p12, p2 = 1, entry = entry1)
  } else if (transition$family == "piecewise exponential") {
    getPCWDistr(U = U1, haz = h12, pw = pw12, t_0 = entry1)
  }
  exit1 <- entry1 + wait_time1

  # Add censoring.
  censTime1 <- simDataOne$censTime
  to1 <- ifelse(censTime1 < exit1, "cens", to1)
  exit1 <- pmin(censTime1, exit1)

  data.frame(
    id = id1,
    from = from1,
    to = to1,
    entry = entry1,
    exit = exit1,
    censTime = censTime1,
    stringsAsFactors = FALSE
  )
}

#' Staggered Study Entry
#'
#' This function adds staggered study entry times to a simulated data set with illness-death model structure.
#'
#' @param simData (`data.frame`)\cr simulated data frame containing entry and exit times
#'   at individual study time scale. See [getSimulatedData()] for details.
#' @param N (`int`)\cr number of patients.
#' @param accrualParam (`string`)\cr possible values are 'time' or 'intensity'.
#' @param accrualValue  (`number`)\cr specifies the accrual intensity. For `accrualParam` equal time,
#'   it describes the length of the accrual period. For `accrualParam` equal intensity, it describes
#'   the number of patients recruited per time unit.  If `accrualValue` is equal to 0,
#'   all patients start at calendar time 0
#'   in the initial state.
#'
#' @return This returns a data set containing a single simulated study containing accrual times,
#'  i.e. staggered study entry.
#'  This is a helper function of [getSimulatedData()].
#'
#' @export
#'
#' @examples
#' simData <- data.frame(
#'   id = c(1, 1, 2, 3), from = c(0, 1, 0, 0), to = c(1, 2, "cens", 2),
#'   entry = c(0, 3, 0, 0),
#'   exit = c(3, 5.3, 5.6, 7.2), censTime = c(6.8, 6.8, 5.6, 9.4)
#' )
#' addStaggeredEntry(simData, 3, accrualParam = "time", accrualValue = 5)
addStaggeredEntry <- function(simData, N, accrualParam, accrualValue) {
  assert_choice(accrualParam, c("time", "intensity"))
  assert_number(accrualValue, lower = 0)

  # Get accrual times in calendar time per individual.
  # If no staggered study entry is present, all individuals have the same entry time 0.
  entry_act <- if (accrualValue != 0) {
    if (accrualParam == "time") {
      stats::runif(N, 0, accrualValue)
    } else if (accrualParam == "intensity") {
      accrualTime <- N / accrualValue
      stats::runif(N, 0, accrualTime)
    }
  } else if (accrualValue == 0) {
    rep(0, N)
  }
  entryAct <- cbind(id = 1:N, entry_act)
  # Combine simulated data with actual entry time.
  # Generate actual times  (actual entry time + individual time).
  simData <- merge(simData, entryAct)
  simData$entryAct <- simData$entry + simData$entry_act
  simData$exitAct <- simData$exit + simData$entry_act
  simData$censAct <- simData$censTime + simData$entry_act
  # Delete temporary helper columns.
  simData$entry_act <- NULL
  simData[, c("censTime")] <- list(NULL)
  simData
}

#' Simulate Data Set from an Illness-Death Model
#'
#' This function creates a single simulated data set for a single treatment arm. It simulates data
#'  from an illness-death model with one row per transition and subject.
#'
#' @param N (`int`)\cr number of patients.
#' @param transition (`TransitionParameters`)\cr transition parameters comprising
#'   `hazards`, corresponding `intervals` and `weibull_rates`, see [exponential_transition()], [piecewise_exponential()]
#'   and [weibull_transition()] for details.
#' @param dropout (`list`)\cr specifies drop-out probability.
#'   Random censoring times are generated using exponential distribution. `dropout$rate` specifies
#'   the drop-out probability per `dropout$time` time units.
#'   If `dropout$rate` is equal to 0, then no censoring is applied.
#' @param accrual (`list`)\cr specifies accrual intensity. See [addStaggeredEntry()] for details.
#'
#' @return This returns a data frame with one row per transition per individual.
#' @details The output data set contains the following columns:
#' - id (`integer`): patient id.
#' - from (`numeric`): starting state of the transition.
#' - to (`character`): final state of the transition.
#' - entry (`numeric`): entry time of the transition on the individual time scale.
#' - exit (`numeric`): exit time of the transition on the individual time scale.
#' - entryAct (`numeric`): entry time of the transition on study time scale.
#' - exitAct (`numeric`):  exit time  of the transition on study time scale.
#' - censAct (`numeric`):  censoring time of the individual on study time scale.
#' @export
#'
#' @examples
#' getSimulatedData(
#'   N = 10,
#'   transition = exponential_transition(h01 = 1, h02 = 1.5, h12 = 1),
#'   dropout = list(rate = 0.3, time = 1),
#'   accrual = list(param = "time", value = 5)
#' )
getSimulatedData <- function(N,
                             transition = exponential_transition(h01 = 1, h02 = 1, h12 = 1),
                             dropout = list(rate = 0, time = 12),
                             accrual = list(param = "time", value = 0)) {
  # Check input parameters.
  assert_int(N, lower = 1L)
  assert_class(transition, "TransitionParameters")
  assert_list(dropout)

  assert(check_number(dropout$rate, lower = 0, upper = 1),
    check_number(dropout$time),
    check_true(dropout$time > 0, ),
    combine = "and", .var.name = "dropout"
  )
  assert_list(accrual)
  # Initialize transition hazards, Weibull rates and intervals.
  h <- transition$hazard
  p <- transition$weibull_rates
  pw <- transition$intervals

  # Get rate parameter for exponential distributed censoring times.
  # if rate=0, no censoring is applied.
  censRate <- if (dropout$rate > 0) {
    -log(1 - dropout$rate) / dropout$time
  } else {
    0
  }
  # Censoring times for all individuals, infinity if no censoring is applied.
  censTime <- if (dropout$rate > 0) {
    stats::rexp(N, censRate)
  } else {
    rep(Inf, N)
  }

  # All individuals start in the initial state.
  entry <- rep(0, N)
  from <- rep(0, N)
  # Waiting time in the initial state 0.
  U <- stats::runif(N)
  # Exponential or Weibull distributed survival times.
  if (transition$family %in% c("exponential", "Weibull")) {
    # For distribution of waiting time in the initial state the all-cause hazard (h01+h02) is needed.
    wait_time <- getWaitTimeSum(U, h$h01, h$h02, p$p01, p$p02, entry)
    # A binomial experiment decides on death or progression.
    numerator <- p$p01 * h$h01 * wait_time^(p$p01 - 1) # nolint
    denumerator <- numerator + p$p02 * h$h02 * wait_time^(p$p02 - 1) # nolint

    # Piecewise exponential distributed survival times.
  } else if (transition$family == "piecewise exponential") {
    Sum_PCW <- getSumPCW(h$h01, h$h02, pw$pw01, pw$pw02)
    wait_time <- getPCWDistr(U, Sum_PCW$hazards, Sum_PCW$intervals, entry)
    # A binomial experiment decides on death or progression.
    numerator <- getPWCHazard(h$h01, pw$pw01, wait_time)
    denumerator <- getPWCHazard(h$h01, pw$pw01, wait_time) +
      getPWCHazard(h$h02, pw$pw02, wait_time)
  }

  to_prob <- stats::rbinom(N, 1, numerator / denumerator)
  to <- ifelse(to_prob == 0, 2, 1)
  exit <- wait_time

  # Add censoring.
  to <- ifelse(censTime < wait_time, "cens", to)
  exit <- pmin(censTime, wait_time)

  simData <- data.frame(
    id = 1:N, from = from, to = to, entry = entry, exit = exit,
    censTime = censTime, stringsAsFactors = FALSE
  )
  is_intermediate_state <- simData$to == 1
  if (any(is_intermediate_state)) {
    # Add 1 -> 2 transition only for patients that are in the intermediate state 1.
    simDataOne <- simData[is_intermediate_state, ]
    newrows <- getOneToTwoRows(simDataOne, transition)
    simData <- rbind(simData, newrows)
    simData <- simData[order(simData$id), ]
  }
  # Add staggered study entry, i.e. study entry at calendar time.
  simData <- addStaggeredEntry(
    simData = simData,
    N = N,
    accrualParam = accrual$param,
    accrualValue = accrual$value
  )
  simData$to <- as.character(simData$to)
  simData
}
