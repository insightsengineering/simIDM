# prepareData ----

test_that("prepareData works as expected", {
  # Create simIDM data for the 4 possible different transition scenarios.
  colnames <- c(
    "id", "trt", "PFStime", "CensoredPFS", "PFSevent", "OStime",
    "CensoredOS", "OSevent", "recruitTime", "OStimeCal", "PFStimeCal"
  )
  patCens1 <- c(1, 1, 0.2, 1, 0, 0.2, 1, 0, 0.1, 0.3, 0.3)
  patCens2 <- c(2, 1, 0.7, 0, 1, 0.8, 1, 0, 1.2, 2, 1.9)
  pat13 <- c(3, 1, 0.4, 0, 1, 0.4, 0, 1, 0.4, 0.8, 0.8)
  pat123 <- c(4, 1, 0.1, 0, 1, 0.25, 0, 1, 0, 0.25, 0.1)
  df <- setNames(data.frame(rbind(patCens1, patCens2, pat13, pat123)), nm = colnames)
  actual <- prepareData(df)
  expect_snapshot(actual)
})

# negLogLik ----

test_that("negLogLik works as expected for Exponential", {
  transition <- exponential_transition(2, 1.3, 0.8)
  data <- prepareData(getClinicalTrials(
    nRep = 1, nPat = 50, seed = 1234, datType = "1rowPatient",
    transitionByArm = list(transition), dropout = list(rate = 0.3, time = 12),
    accrual = list(param = "intensity", value = 7)
  )[[1]])
  actual1 <- negLogLik(transition, data)
  expect_equal(actual1, 58.65772)
})

test_that("negLogLik works as expected for Weibull", {
  transition <- weibull_transition(h01 = 0.2, h02 = 0.5, h12 = 1.6, p01 = 1, p02 = 2.5, p12 = 3)
  data <- prepareData(getClinicalTrials(
    nRep = 1, nPat = 50, seed = 1234, datType = "1rowPatient",
    transitionByArm = list(transition), dropout = list(rate = 0.3, time = 12),
    accrual = list(param = "intensity", value = 7)
  )[[1]])
  actual2 <- negLogLik(transition, data)
  expect_equal(actual2, 54.644006)
})

# haz ----

test_that("haz works as expected for Exponential", {
  transition <- exponential_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6)
  actual1 <- haz(transition, 0.4, 2)
  expect_equal(actual1, 1.5)
})

test_that("haz works as expected for Weibull", {
  transition <- weibull_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6, p01 = 2, p02 = 2.5, p12 = 3)
  actual2 <- haz(transition, 0.4, 2)
  expect_equal(actual2, 0.9486833)
})

# survTrans----

test_that("survTrans works as expected for Exponential", {
  transition <- exponential_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6)
  actual1 <- survTrans(transition, 0.4, 2)
  expect_equal(actual1, 0.54881164)
})

test_that("survTrans works as expected for Weibull", {
  transition <- weibull_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6, p01 = 2, p02 = 2.5, p12 = 3)
  actual2 <- survTrans(transition, 0.1, 3)
  expect_equal(actual2, 0.99840128)
})

# getInit ----

test_that("getInit works as expected for Exponential", {
  transition <- exponential_transition(h01 = 2.2, h02 = 0.5, h12 = 1.3)
  actual1 <- getInit(transition)
  expect_equal(actual1, c(2.2, 0.5, 1.3))
})

test_that("getInit works as expected for Weibull", {
  transition <- weibull_transition(h01 = 0.2, h02 = 0.5, h12 = 1.6, p01 = 1, p02 = 2.5, p12 = 3)
  actual2 <- getInit(transition)
  expect_equal(actual2, c(0.2, 0.5, 1.6, 1, 2.5, 3))
})

# getTarget ----

test_that("getTarget works as expected for Exponential", {
  transition <- exponential_transition(2, 1.3, 0.8)
  data <- prepareData(getClinicalTrials(
    nRep = 1, nPat = 50, seed = 1234, datType = "1rowPatient",
    transitionByArm = list(transition), dropout = list(rate = 0.3, time = 12),
    accrual = list(param = "intensity", value = 7)
  )[[1]])
  params <- c(1.2, 1.5, 1.6)
  target <- getTarget(transition)
  actual1 <- target(params, data)
  expect_equal(actual1, 84.68301)
})

test_that("getTarget works as expected for Weibull", {
  transition <- weibull_transition(h01 = 0.2, h02 = 0.5, h12 = 1.6, p01 = 1, p02 = 2.5, p12 = 3)
  data <- prepareData(getClinicalTrials(
    nRep = 1, nPat = 50, seed = 1234, datType = "1rowPatient",
    transitionByArm = list(transition), dropout = list(rate = 0.3, time = 12),
    accrual = list(param = "intensity", value = 7)
  )[[1]])
  params <- c(1.2, 1.5, 1.6, 2, 1, 2)
  target <- getTarget(transition)
  actual2 <- target(params, data)
  expect_equal(actual2, 103.6357444)
})

# getResults ----

test_that("getResults works as expected for Exponential", {
  results <- c(1.2, 1.5, 1.6)
  actual1 <- getResults(exponential_transition(), results)
  expect_identical(actual1$hazards, list(h01 = 1.2, h02 = 1.5, h12 = 1.6))
})

test_that("getResults works as expected for weibull", {
  results <- c(1.2, 1.5, 1.6, 2, 1, 0.5)
  actual2 <- getResults(weibull_transition(), results)
  expect_identical(actual2$hazards, list(h01 = 1.2, h02 = 1.5, h12 = 1.6))
  expect_identical(actual2$weibull_rates, list(p01 = 2, p02 = 1, p12 = 0.5))
})

# estimateParams ----

test_that("estimateParams estimates the true parameters correctly for Exponential", {
  transition <- exponential_transition(2, 1.3, 0.8)
  data <- getClinicalTrials(
    nRep = 1, nPat = 100000, seed = 123, datType = "1rowPatient",
    transitionByArm = list(transition), dropout = list(rate = 0.3, time = 1),
    accrual = list(param = "intensity", value = 500)
  )[[1]]
  actual1 <- estimateParams(data, transition)
  expect_equal(actual1$hazards, list(h01 = 2, h02 = 1.3, h12 = 0.8), tolerance = 1e-2)
  expect_identical(actual1$weibull_rates, list(p01 = 1, p02 = 1, p12 = 1))
  expect_identical(class(actual1), c("ExponentialTransition", "TransitionParameters"))
})

test_that("estimateParams estimates the true parameters correctly for Weibull", {
  transition <- weibull_transition(h01 = 0.4, h02 = 0.9, h12 = 1.6, p01 = 1, p02 = 0.5, p12 = 1.9)
  data <- getClinicalTrials(
    nRep = 1, nPat = 100000, seed = 123, datType = "1rowPatient",
    transitionByArm = list(transition), dropout = list(rate = 0.3, time = 1),
    accrual = list(param = "intensity", value = 500)
  )[[1]]
  actual2 <- estimateParams(data, transition)
  expect_equal(actual2$hazards, list(h01 = 0.4, h02 = 0.9, h12 = 1.6), tolerance = 1e-2)
  expect_equal(actual2$weibull_rates, list(p01 = 1, p02 = 0.5, p12 = 1.9), tolerance = 1e-2)
  expect_identical(class(actual2), c("WeibullTransition", "TransitionParameters"))
})
