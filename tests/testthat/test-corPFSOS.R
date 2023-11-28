# corPFSOS ----

test_that("corPFSOS returns correct central estimate", {
  transition <- exponential_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6)
  data <- getClinicalTrials(
    nRep = 1, nPat = c(20000), seed = 1234, datType = "1rowTransition",
    transitionByArm = list(transition), dropout = list(rate = 0.5, time = 12),
    accrual = list(param = "intensity", value = 7)
  )[[1]]
  actual <- corPFSOS(data, transition = exponential_transition(), bootstrap = FALSE)
  expect_equal(actual[[1]], corTrans(transition = transition), tolerance = 1e-2)
})
