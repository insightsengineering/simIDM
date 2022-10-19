# getSimulatedData --- ----
test_that("getSimulatedData generates distributions as expected", {
  transition1 <- exponential_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6)
  transition2 <- exponential_transition(h01 = 1, h02 = 1.3, h12 = 1.7)
  actual <- getClinicalTrials(
    nRep = 1000, nPat = c(200, 200), seed = 1234, datType = "1rowTransition",
    transitionByArm = list(transition1, transition2), dropout = list(rate = 0, time = 12),
    accrual = list(param = "intensity", value = 7)
  )
  # commpare true Nelson-Aalen estimator with simulated one.

  # transition matrix
  tra <- matrix(ncol = 3, nrow = 3, FALSE)
  tra[1, 2:3] <- TRUE
  tra[2, c(3)] <- TRUE

  # Nelson-Aalen estimator
  times <- seq(0, 3, 0.01)
  library(mvna)
  EstimatedNA <- lapply(seq_along(actual), function(j) {
    na <- NULL
    na <- mvna(actual[[j]][actual[[j]]$trt == 1, ], c("0", "1", "2"), tra, "cens")
    na_predict <- predict(na, times,
      tr.choice = c("0 1"),
      level = 0.95, var.type = c("aalen"), ci.fun = c("log")
    )[["0 1"]][, "na"]
    return(na_predict)
  })
  EstimatedNAMean <- rowMeans(matrix(unlist(EstimatedNA), ncol = length(EstimatedNA)))

  # true NA - trt 1 0 -> 1 transition

  trueNA <- transition1$hazards$h01 * times

  if (interactive()) {
    plot(times, EstimatedNAMean, type = "l")
    lines(times, trueNA, col = "red")
  }
  tol1 <- 0.01
  expect_true(all(abs(EstimatedNAMean[1:130] - trueNA[1:130]) <= tol1))
})
