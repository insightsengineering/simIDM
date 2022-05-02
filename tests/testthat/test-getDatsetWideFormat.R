# getDatasetWideFormat  ----

test_that("getDatasetWideFormat works as expected", {
  transition1 <- weibull_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6, p01 = 0.8, p02 = 0.9, p12 = 1)
  transition2 <- weibull_transition(h01 = 1, h02 = 1.3, h12 = 1.7, p01 = 1.1, p02 = 0.9, p12 = 1.1)

  simStudy <- getOneClinicalTrial(
    nPat = c(30, 30), transitionByArm = list(transition1, transition2),
    dropout = list(rate = 0.5, time = 12),
    accrual = list(param = "time", value = 0)
  )

  actual <- getDatasetWideFormat(simStudy)
  expect_true(length(actual$id[actual$trt == 1]) == 30)
  expect_true(length(actual$id[actual$trt == 2]) == 30)

  expect_true(all(actual$PFSevent + actual$CensoredPFS == 1))
  expect_true(all(actual$OSevent + actual$CensoredOS == 1))

  expect_true(all(actual$OStime <= actual$OStimeCal))
  expect_true(all(actual$PFStime <= actual$PFStimeCal))

  expect_named(actual, c(
    "id", "trt", "PFStime", "CensoredPFS", "PFSevent", "OStime", "CensoredOS",
    "OSevent", "recruitTime", "OStimeCal", "PFStimeCal"
  ))
})
