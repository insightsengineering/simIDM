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

# powerEmp ----

test_that("powerEmp works as expected", {
  transition1 <- exponential_transition(h01 = 0.06, h02 = 0.3, h12 = 0.3)
  transition2 <- exponential_transition(h01 = 0.1, h02 = 0.4, h12 = 0.3)
  simTrials <- getClinicalTrials(
    nRep = 50, nPat = c(800, 800), seed = 1234, datType = "1rowPatient",
    transitionByArm = list(transition1, transition2), dropout = list(rate = 0.5, time = 12),
    accrual = list(param = "intensity", value = 7)
  )
  actual <- powerEmp(
    simTrials = simTrials, criticalPFS = 2.4, criticalOS = 2.2,
    eventNumPFS = 300, eventNumOS = 500
  )
  expect_equal(actual, c(
    "powerPFS" = 0.74,
    "powerOS" = 0.52,
    "powerAtLeast" = 0.78,
    "powerJoint" = 0.48
  ))
})
