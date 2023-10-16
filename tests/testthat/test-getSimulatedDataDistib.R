# helper function for tests.
library(mvna)
getEstimatedNA <- function(actual, times) {
  # transition matrix
  tra <- matrix(ncol = 3, nrow = 3, FALSE)
  tra[1, 2:3] <- TRUE
  tra[2, c(3)] <- TRUE

  # Nelson-Aalen estimator.
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
  return(EstimatedNAMean)
}



# getSimulatedData --- ----
test_that("getSimulatedData generates distributions as expected - Exponential", {
  transition1 <- exponential_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6)
  transition2 <- exponential_transition(h01 = 1, h02 = 1.3, h12 = 1.7)
  actual <- getClinicalTrials(
    nRep = 1000, nPat = c(200, 200), seed = 1234, datType = "1rowTransition",
    transitionByArm = list(transition1, transition2), dropout = list(rate = 0, time = 12),
    accrual = list(param = "intensity", value = 7)
  )
  # compare true Nelson-Aalen estimator with simulated one.
  times <- seq(0, 3, 0.01)
  EstimatedNAMean <- getEstimatedNA(actual, times)

  # true NA - trt 1 0 -> 1 transition.

  trueNA <- transition1$hazards$h01 * times

  if (interactive()) {
    plot(times, EstimatedNAMean, type = "l")
    lines(times, trueNA, col = "red")
  }
  tol1 <- 0.01
  expect_true(all(abs(EstimatedNAMean[1:130] - trueNA[1:130]) <= tol1))
})




# getSimulatedData --- ----
test_that("getSimulatedData generates distributions as expected - Weibull", {
  transition1 <- weibull_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6, p01 = 1.2, p02 = 1, p12 = 0.3)
  transition2 <- weibull_transition(h01 = 1, h02 = 1.3, h12 = 1.7, p01 = 1.2, p02 = 1, p12 = 0.3)
  actual <- getClinicalTrials(
    nRep = 1000, nPat = c(200, 200), seed = 1234, datType = "1rowTransition",
    transitionByArm = list(transition1, transition2), dropout = list(rate = 0, time = 12),
    accrual = list(param = "intensity", value = 7)
  )
  # compare true Nelson-Aalen estimator with simulated one.
  times <- seq(0, 3, 0.01)
  EstimatedNAMean <- getEstimatedNA(actual, times)

  # true NA - trt 1 0 -> 1 transition.

  trueNA <- transition1$hazards$h01 * times^transition1$weibull_rates$p01

  if (interactive()) {
    plot(times, EstimatedNAMean, type = "l")
    lines(times, trueNA, col = "red")
  }
  tol1 <- 0.01
  expect_true(all(abs(EstimatedNAMean[1:100] - trueNA[1:100]) <= tol1))
})





# getSimulatedData --- ----
test_that("getSimulatedData generates distributions as expected - PW", {
  transition1 <- piecewise_exponential(
    h01 = c(1, 0.8, 1.3), h02 = c(1.6, 1.5, 1), h12 = c(1, 1.6, 1),
    pw01 = c(0, 0.5, 1), pw02 = c(0, 3, 8), pw12 = c(0, 2, 4)
  )
  transition2 <- piecewise_exponential(
    h01 = c(1, 0.9, 1.5), h02 = c(1.7, 1.8, 1), h12 = c(1.2, 1.6, 1),
    pw01 = c(0, 1, 7), pw02 = c(0, 3, 8), pw12 = c(0, 2, 4)
  )
  actual <- getClinicalTrials(
    nRep = 1000, nPat = c(200, 200), seed = 1234, datType = "1rowTransition",
    transitionByArm = list(transition1, transition2), dropout = list(rate = 0, time = 12),
    accrual = list(param = "intensity", value = 7)
  )
  # compare true Nelson-Aalen estimator with simulated one.
  times <- seq(0, 3, 0.01)
  EstimatedNAMean <- getEstimatedNA(actual, times)

  # true NA - trt 1 0 -> 1 transition.
  pw01TRT1 <- transition1$intervals$pw01
  h01TRT1 <- transition1$hazards$h01
  trueNA <- truePWC <- sapply(times, function(t) {
    if (t <= pw01TRT1[2]) {
      return(h01TRT1[1] * (t - pw01TRT1[1]))
    } else if (t <= pw01TRT1[3]) {
      return(h01TRT1[1] * (pw01TRT1[2] - pw01TRT1[1]) + h01TRT1[2] * (t - pw01TRT1[2]))
    } else {
      return(h01TRT1[1] * (pw01TRT1[2] - pw01TRT1[1]) + h01TRT1[2] *
        (pw01TRT1[3] - pw01TRT1[2]) + h01TRT1[3] * (t - pw01TRT1[3]))
    }
  })


  if (interactive()) {
    plot(times, EstimatedNAMean, type = "l")
    lines(times, trueNA, col = "red")
  }
  tol1 <- 0.01
  expect_true(all(abs(EstimatedNAMean[1:100] - trueNA[1:100]) <= tol1))
})
