# getOneClinicalTrial ----

test_that("getOneClinicalTrial works as expected", {
  transition1 <- weibull_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6, p01 = 0.8, p02 = 0.9, p12 = 1)
  transition2 <- weibull_transition(h01 = 1, h02 = 1.3, h12 = 1.7, p01 = 1.1, p02 = 0.9, p12 = 1.1)
  transition3 <- weibull_transition(h01 = 1.1, h02 = 1, h12 = 1.5, p01 = 0.9, p02 = 0.9, p12 = 1.2)

  actual <- getOneClinicalTrial(
    nPat = c(30, 20, 40), transitionByArm = list(transition1, transition2, transition3),
    dropout = list(rate = 0.1, time = 12),
    accrual = list(param = "intensity", value = 10)
  )
  actualUnique <- actual[!duplicated(actual$id), ]
  expect_true(length(actualUnique$id[actualUnique$trt == 1]) == 30)
  expect_true(length(actualUnique$id[actualUnique$trt == 2]) == 20)
  expect_true(length(actualUnique$id[actualUnique$trt == 3]) == 40)
})


test_that("getOneClinicalTrial works as expected - accrual and dropout as list", {
  transition1 <- weibull_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6, p01 = 0.8, p02 = 0.9, p12 = 1)
  transition2 <- weibull_transition(h01 = 1, h02 = 1.3, h12 = 1.7, p01 = 1.1, p02 = 0.9, p12 = 1.1)
  transition3 <- weibull_transition(h01 = 1.1, h02 = 1, h12 = 1.5, p01 = 0.9, p02 = 0.9, p12 = 1.2)

  dropout1 <- list(rate = 0.1, time = 12)
  dropout2 <- list(rate = 0.2, time = 12)
  dropout3 <- list(rate = 0.3, time = 12)
  accrual1 <- list(param = "intensity", value = 10)
  accrual2 <- list(param = "intensity", value = 5)
  accrual3 <- list(param = "intensity", value = 15)

  actual <- getOneClinicalTrial(
    nPat = c(30, 20, 40), transitionByArm = list(transition1, transition2, transition3),
    dropout = list(dropout1, dropout2, dropout3),
    accrual = list(accrual1, accrual2, accrual3)
  )
  actualUnique <- actual[!duplicated(actual$id), ]
  expect_true(length(actualUnique$id[actualUnique$trt == 1]) == 30)
  expect_true(length(actualUnique$id[actualUnique$trt == 2]) == 20)
  expect_true(length(actualUnique$id[actualUnique$trt == 3]) == 40)
})

test_that("getOneClinicalTrial creates expected data set", {
  set.seed(1243)
  transition1 <- weibull_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6, p01 = 0.8, p02 = 0.9, p12 = 1)
  transition2 <- weibull_transition(h01 = 1, h02 = 1.3, h12 = 1.7, p01 = 1.1, p02 = 0.9, p12 = 1.1)

  actual <- getOneClinicalTrial(
    nPat = c(30, 30), transitionByArm = list(transition1, transition2),
    dropout = list(rate = 0.5, time = 12),
    accrual = list(param = "time", value = 0)
  )
  row1 <- data.frame(
    id = 1, from = 0, to = "1", entry = 0.0000000, exit = 0.186746130343395766,
    entryAct = 0, exitAct = 0.186746130343395766, censAct = 19.10280613143678, trt = 1,
    stringsAsFactors = FALSE
  )
  expect_equal(actual[1, ], row1)
})

test_that("getOneClinicalTrial also works in edge case", {
  n_pat <- c(1, 1)
  param_control <- c(0.04781994, 0.03882345, 0.31313900)
  param_i <- c(0.06971534, 0.00000000, 0.24160661)
  sep_time <- 2
  transition_i <- piecewise_exponential(
    pw01 = c(0, sep_time),
    pw02 = c(0, sep_time),
    pw12 = c(0, sep_time),
    h01 = c(param_control[1], param_i[1]),
    h02 = c(param_control[2], param_i[2]),
    h12 = c(param_control[3], param_i[3])
  )
  transition_control <- piecewise_exponential(
    pw01 = c(0, sep_time),
    pw02 = c(0, sep_time),
    pw12 = c(0, sep_time),
    h01 = rep(param_control[1], 2),
    h02 = rep(param_control[2], 2),
    h12 = rep(param_control[3], 2)
  )
  transition_by_arm <- list(transition_i, transition_control)
  dropout <- list(rate = 0, time = 1)
  accrual <- list(param = "intensity", value = 5)
  set.seed(977)
  result <- expect_silent(getOneClinicalTrial(
    nPat = n_pat,
    transitionByArm = transition_by_arm,
    dropout = dropout,
    accrual = accrual
  ))
  expect_true(any(result$to == 1))
  expect_true(any(result$from == 1))
})

test_that("getOneClinicalTrial gives a warning if there are no progression to death transitions at all", {
  n_pat <- c(1, 1)
  param_control <- c(0.04781994, 0.03882345, 0.31313900)
  param_i <- c(0.06971534, 0.00000000, 0.24160661)
  sep_time <- 2
  transition_i <- piecewise_exponential(
    pw01 = c(0, sep_time),
    pw02 = c(0, sep_time),
    pw12 = c(0, sep_time),
    h01 = c(param_control[1], param_i[1]),
    h02 = c(param_control[2], param_i[2]),
    h12 = c(param_control[3], param_i[3])
  )
  transition_control <- piecewise_exponential(
    pw01 = c(0, sep_time),
    pw02 = c(0, sep_time),
    pw12 = c(0, sep_time),
    h01 = rep(param_control[1], 2),
    h02 = rep(param_control[2], 2),
    h12 = rep(param_control[3], 2)
  )
  transition_by_arm <- list(transition_i, transition_control)
  dropout <- list(rate = 0, time = 1)
  accrual <- list(param = "intensity", value = 5)
  set.seed(59)
  result <- expect_warning(
    getOneClinicalTrial(
      nPat = n_pat,
      transitionByArm = transition_by_arm,
      dropout = dropout,
      accrual = accrual
    ),
    "no progression to death transitions included in the simulated data"
  )
  expect_false(any(result$to == 1))
  expect_false(any(result$from == 1))
})

# getClinicalTrials ----

test_that("getClinicalTrials works as expected if it creates 1 row per patient", {
  transition1 <- exponential_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6)
  transition2 <- exponential_transition(h01 = 1, h02 = 1.3, h12 = 1.7)
  actual <- getClinicalTrials(
    nRep = 10, nPat = c(20, 20), seed = 1234, datType = "1rowPatient",
    transitionByArm = list(transition1, transition2), dropout = list(rate = 0.5, time = 12),
    accrual = list(param = "intensity", value = 7)
  )
  expect_type(actual, "list")
  expect_length(actual, 10)
  # 1 row per patients?
  expect_true(all(sapply(actual, nrow) == 40))
})

test_that("getClinicalTrials works as expected if accrual and dropout are lists", {
  transition1 <- exponential_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6)
  transition2 <- exponential_transition(h01 = 1, h02 = 1.3, h12 = 1.7)
  dropout1 <- list(rate = 0.1, time = 12)
  dropout2 <- list(rate = 0.2, time = 12)

  accrual1 <- list(param = "intensity", value = 10)
  accrual2 <- list(param = "intensity", value = 5)

  actual <- getClinicalTrials(
    nRep = 10, nPat = c(20, 20), seed = 1234, datType = "1rowPatient",
    transitionByArm = list(transition1, transition2), dropout = list(dropout1, dropout2),
    accrual = list(accrual1, accrual2)
  )
  expect_type(actual, "list")
  expect_length(actual, 10)
  # 1 row per patients?
  expect_true(all(sapply(actual, nrow) == 40))
})

test_that("getClinicalTrials works as expected if it creates 1 row per transition", {
  transition1 <- exponential_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6)
  transition2 <- exponential_transition(h01 = 1, h02 = 1.3, h12 = 1.7)
  actual <- getClinicalTrials(
    nRep = 10, nPat = c(20, 20), seed = 1234, datType = "1rowTransition",
    transitionByArm = list(transition1, transition2), dropout = list(rate = 0.5, time = 12),
    accrual = list(param = "intensity", value = 7)
  )

  expect_type(actual, "list")
  expect_length(actual, 10)
  # at least one row per patient and not more than two per patient?
  expect_true(all(sapply(actual, nrow) >= 40 & sapply(actual, nrow) <= 80))
})

test_that("getClinicalTrials works as expected and generates data according to Exponential", {
  transition <- exponential_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6)
  actual <- getClinicalTrials(
    nRep = 1, nPat = c(2000, 2000), seed = 1234, datType = "1rowPatient",
    transitionByArm = list(transition, transition), dropout = list(rate = 0, time = 12),
    accrual = list(param = "intensity", value = 0)
  )[[1]]
  # Theoretical CDF PFS:
  ExpCDFPFS <- function(t, h01, h02) {
    1 - ExpSurvPFS(t, h01, h02)
  }
  # Theoretical CDF OS:
  ExpCDFOS <- function(t, h01, h02, h12) {
    1 - ExpSurvOS(t, h01, h02, h12)
  }
  # Kolmogorov-Smirnov Tests:
  testPFS <- ks.test(x = actual$PFStime, y = "ExpCDFPFS", h01 = 1.2, h02 = 1.5)$p.value
  testOS <- ks.test(x = actual$OStime, y = "ExpCDFOS", h01 = 1.2, h02 = 1.5, h12 = 1.6)$p.value
  expect_true(testPFS > 0.05)
  expect_true(testOS > 0.05)
})

test_that("getClinicalTrials works as expected and generates data according to Weibull", {
  transition <- weibull_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6, p01 = 2, p02 = 0.4, p12 = 1.3)
  actual <- getClinicalTrials(
    nRep = 1, nPat = c(2000, 2000), seed = 1234, datType = "1rowPatient",
    transitionByArm = list(transition, transition), dropout = list(rate = 0, time = 12),
    accrual = list(param = "intensity", value = 0)
  )[[1]]
  # Theoretical CDF PFS:
  WeibCDFPFS <- function(t, h01, h02, p01, p02) {
    1 - WeibSurvPFS(t, h01, h02, p01, p02)
  }
  # Theoretical CDF OS:
  WeibCDFOS <- function(t, h01, h02, h12, p01, p02, p12) {
    1 - WeibSurvOS(t, h01, h02, h12, p01, p02, p12)
  }
  # Kolmogorov-Smirnov Tests:
  testPFS <- ks.test(
    x = actual$PFStime, y = "WeibCDFPFS", h01 = 1.2, h02 = 1.5,
    p01 = 2, p02 = 0.4
  )$p.value
  testOS <- ks.test(
    x = actual$OStime, y = "WeibCDFOS", h01 = 1.2, h02 = 1.5, h12 = 1.6,
    p01 = 2, p02 = 0.4, p12 = 1.3
  )$p.value
  expect_true(testPFS > 0.05)
  expect_true(testOS > 0.05)
})

test_that("getClinicalTrials works as expected and generates data according to PWC", {
  transition <- piecewise_exponential(
    h01 = c(1, 1.2, 1), h02 = c(0.3, 1, 2), h12 = c(2, 1, 2),
    pw01 = c(0, 3, 8), pw02 = c(0, 6, 7), pw12 = c(0, 8, 9)
  )
  actual <- getClinicalTrials(
    nRep = 1, nPat = c(2000, 2000), seed = 1234, datType = "1rowPatient",
    transitionByArm = list(transition, transition), dropout = list(rate = 0, time = 12),
    accrual = list(param = "intensity", value = 0)
  )[[1]]
  # Theoretical CDF PFS:
  PWCCDFPFS <- function(t, h01, h02, pw01, pw02) {
    1 - PWCsurvPFS(t, h01, h02, pw01, pw02)
  }
  # Theoretical CDF OS:
  PWCCDFOS <- function(t, h01, h02, h12, pw01, pw02, pw12) {
    1 - PWCsurvOS(t, h01, h02, h12, pw01, pw02, pw12)
  }
  # Kolmogorov-Smirnov Tests:
  testPFS <- ks.test(
    x = actual$PFStime, y = "PWCCDFPFS", h01 = c(1, 1.2, 1), h02 = c(0.3, 1, 2),
    pw01 = c(0, 3, 8), pw02 = c(0, 6, 7)
  )$p.value
  testOS <- ks.test(
    x = actual$OStime, y = "PWCCDFOS", h01 = c(1, 1.2, 1), h02 = c(0.3, 1, 2), h12 = c(2, 1, 2),
    pw01 = c(0, 3, 8), pw02 = c(0, 6, 7), pw12 = c(0, 8, 9)
  )$p.value
  expect_true(testPFS > 0.05)
  expect_true(testOS > 0.05)
})
