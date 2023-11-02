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
#' @keywords internal
prepareData <- function(data) {
  assert_data_frame(data, min.cols = 9, max.cols = 11)
  if (ncol(data) == 9) {
    data <- getDatasetWideFormat(data)
  } else {
    assert_subset(c(
      "id", "trt", "PFStime", "CensoredPFS", "PFSevent", "OStime",
      "CensoredOS", "OSevent", "recruitTime",
      "OStimeCal", "PFStimeCal"
    ), names(data))
  }

  # Transform simIDM trial data to log-likelihood-compatible format.
  dataNew <- suppressWarnings(mstate::msprep(
    time = c("recruitTime", "PFStime", "OStime"),
    status = c("trt", "PFSevent", "OSevent"),
    data = data,
    trans = mstate::trans.illdeath(),
    id = data$id
  ))
  names(dataNew)[5:6] <- c("entry", "exit")
  # Correct msprep results for uncensored PFS=OS events.
  ids <- data$id[data$PFStimeCal == data$OStimeCal & data$CensoredPFS == 0]
  dataNew <- dataNew[!(dataNew$id %in% ids & dataNew$trans == 3), ]
  dataNew$status[dataNew$id %in% ids] <- abs(dataNew$status[dataNew$id %in% ids] - 1)

  as.data.frame(dataNew[, -7], row.names = seq_len(nrow(dataNew)))
}
