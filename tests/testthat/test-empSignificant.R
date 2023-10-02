# logRankTest ----

test_that("logRankTest works as expected", {
  transition1 <- exponential_transition(h01 = 0.06, h02 = 0.3, h12 = 0.3)
  transition2 <- exponential_transition(h01 = 0.1, h02 = 0.4, h12 = 0.3)
  simTrial <- getClinicalTrials(
    nRep = 1, nPat = c(800, 800), seed = 1234, datType = "1rowPatient",
    transitionByArm = list(transition1, transition2), dropout = list(rate = 0.5, time = 12),
    accrual = list(param = "intensity", value = 7)
  )[[1]]
  actual <- logRankTest(data = simTrial, typeEvent = "OS", critical = 3.4)
  expect_equal(actual, TRUE)

  actual2 <- logRankTest(data = simTrial, typeEvent = "PFS", critical = 6)
  expect_equal(actual2, FALSE)
})

# passedLogRank ----

test_that("passedLogRank works as expected", {
  transition1 <- exponential_transition(h01 = 0.06, h02 = 0.3, h12 = 0.3)
  transition2 <- exponential_transition(h01 = 0.1, h02 = 0.4, h12 = 0.3)
  simTrials <- getClinicalTrials(
    nRep = 3, nPat = c(800, 800), seed = 1234, datType = "1rowPatient",
    transitionByArm = list(transition1, transition2), dropout = list(rate = 0.5, time = 12),
    accrual = list(param = "intensity", value = 7)
  )
  actual <- passedLogRank(simTrials = simTrials, typeEvent = "PFS", eventNum = 300, critical = 2.4)
  expect_equal(actual, c(TRUE, TRUE, FALSE))

  actual2 <- passedLogRank(simTrials = simTrials, typeEvent = "OS", eventNum = 300, critical = 2.4)
  expect_equal(actual2, c(FALSE, FALSE, FALSE))
})

# empSignificant ----

test_that("empSignificant works as expected", {
  transition1 <- exponential_transition(h01 = 0.06, h02 = 0.3, h12 = 0.3)
  transition2 <- exponential_transition(h01 = 0.1, h02 = 0.4, h12 = 0.3)
  simTrials <- getClinicalTrials(
    nRep = 50, nPat = c(800, 800), seed = 1234, datType = "1rowPatient",
    transitionByArm = list(transition1, transition2), dropout = list(rate = 0.5, time = 12),
    accrual = list(param = "intensity", value = 7)
  )
  actual <- empSignificant(
    simTrials = simTrials, criticalPFS = 2.4, criticalOS = 2.2,
    eventNumPFS = 300, eventNumOS = 500
  )
  expect_equal(actual, list(
    "significantPFS" = 0.74,
    "significantOS" = 0.52,
    "significantAtLeastOne" = 0.78,
    "significantBoth" = 0.48
  ))
})
