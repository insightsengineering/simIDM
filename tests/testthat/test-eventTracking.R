# getTimePoint ----

test_that("getTimePoint  works as expected", {
  transition1 <- weibull_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6, p01 = 0.8, p02 = 0.9, p12 = 1)
  transition2 <- weibull_transition(h01 = 1, h02 = 1.3, h12 = 1.7, p01 = 1.1, p02 = 0.9, p12 = 1.1)

  simStudy <- getOneClinicalTrial(
    nPat = c(30, 30), transitionByArm = list(transition1, transition2),
    dropout = list(rate = 0.5, time = 12),
    accrual = list(param = "time", value = 0)
  )
  simStudyWide <- getDatasetWideFormat(simStudy)

  actual <- getTimePoint(data = simStudyWide, eventNum = 10, typeEvent = "OS")
  actual2 <- getTimePoint(data = simStudyWide, eventNum = 6, typeEvent = "PFS")
  expect_equal(length(simStudyWide$id[simStudyWide$OSevent == 1 & simStudyWide$OStime <= actual]), 10)
  expect_equal(length(simStudyWide$id[simStudyWide$PFSevent == 1 & simStudyWide$PFStime <= actual2]), 6)
})


test_that("getTimePoint works as expected by Arm", {
  transition1 <- weibull_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6, p01 = 0.8, p02 = 0.9, p12 = 1)
  transition2 <- weibull_transition(h01 = 1, h02 = 1.3, h12 = 1.7, p01 = 1.1, p02 = 0.9, p12 = 1.1)

  simStudy <- getOneClinicalTrial(
    nPat = c(30, 30), transitionByArm = list(transition1, transition2),
    dropout = list(rate = 0.5, time = 12),
    accrual = list(param = "time", value = 0)
  )
  simStudyWide <- getDatasetWideFormat(simStudy)

  actual <- getTimePoint(data = simStudyWide, eventNum = 5, typeEvent = "OS", byArm = TRUE)
  actual2 <- getTimePoint(data = simStudyWide, eventNum = 7, typeEvent = "PFS", byArm = TRUE)
  expect_equal(length(simStudyWide$id[simStudyWide$OSevent == 1 & simStudyWide$OStime <= actual[1] &
    simStudyWide$trt == 1]), 5)
  expect_equal(length(simStudyWide$id[simStudyWide$OSevent == 1 & simStudyWide$OStime <= actual[2] &
    simStudyWide$trt == 2]), 5)
  expect_equal(length(simStudyWide$id[simStudyWide$PFSevent == 1 & simStudyWide$PFStime <= actual2[1] &
    simStudyWide$trt == 1]), 7)
  expect_equal(length(simStudyWide$id[simStudyWide$PFSevent == 1 & simStudyWide$PFStime <= actual2[2] &
    simStudyWide$trt == 2]), 7)
})


# censoringByNumberEvents----

test_that("censoringByNumberEvents works as expected", {
  transition1 <- weibull_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6, p01 = 0.8, p02 = 0.9, p12 = 1)
  transition2 <- weibull_transition(h01 = 1, h02 = 1.3, h12 = 1.7, p01 = 1.1, p02 = 0.9, p12 = 1.1)

  simStudy <- getOneClinicalTrial(
    nPat = c(30, 30), transitionByArm = list(transition1, transition2),
    dropout = list(rate = 0.5, time = 12),
    accrual = list(param = "time", value = 0)
  )
  simStudyWide <- getDatasetWideFormat(simStudy)

  actual <- censoringByNumberEvents(data = simStudyWide, eventNum = 12, typeEvent = "OS")
  actual2 <- censoringByNumberEvents(data = simStudyWide, eventNum = 16, typeEvent = "PFS")
  expect_equal(length(simStudyWide$id[actual$OSevent == 1]), 12)
  expect_equal(length(simStudyWide$id[actual2$PFSevent == 1]), 16)
})

# trackEventsPerTrial----

test_that("trackEventsPerTrial works as expected by arm", {
  transition1 <- weibull_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6, p01 = 0.8, p02 = 0.9, p12 = 1)
  transition2 <- weibull_transition(h01 = 1, h02 = 1.3, h12 = 1.7, p01 = 1.1, p02 = 0.9, p12 = 1.1)

  simStudy <- getOneClinicalTrial(
    nPat = c(30, 30), transitionByArm = list(transition1, transition2),
    dropout = list(rate = 0.5, time = 2),
    accrual = list(param = "time", value = 0)
  )
  simStudyWide <- getDatasetWideFormat(simStudy)

  actual <- trackEventsPerTrial(data = simStudyWide, timeP = 0.3, byArm = TRUE)
  expect_equal(actual[["1"]]["Recruited", ] - (actual[["1"]]["OS", ] + actual[["1"]]["Censored", ]
    + actual[["1"]]["Ongoing", ]), 0)
  expect_true(actual[["1"]]["PFS", ] >= actual[["1"]]["OS", ])

  expect_equal(actual[["2"]]["Recruited", ] - (actual[["2"]]["OS", ] + actual[["2"]]["Censored", ]
    + actual[["2"]]["Ongoing", ]), 0)
  expect_true(actual[["2"]]["PFS", ] >= actual[["2"]]["OS", ])
})
