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
    accrual = list(param = "time", value = 5)
  )
  simStudyWide <- getDatasetWideFormat(simStudy)

  actual <- getTimePoint(data = simStudyWide, eventNum = 5, typeEvent = "OS", byArm = TRUE)
  actual2 <- getTimePoint(data = simStudyWide, eventNum = 7, typeEvent = "PFS", byArm = TRUE)
  expect_equal(length(simStudyWide$id[simStudyWide$OSevent == 1 & simStudyWide$OStimeCal <= actual[1] &
    simStudyWide$trt == 1]), 5)
  expect_equal(length(simStudyWide$id[simStudyWide$OSevent == 1 & simStudyWide$OStimeCal <= actual[2] &
    simStudyWide$trt == 2]), 5)
  expect_equal(length(simStudyWide$id[simStudyWide$PFSevent == 1 & simStudyWide$PFStimeCal <= actual2[1] &
    simStudyWide$trt == 1]), 7)
  expect_equal(length(simStudyWide$id[simStudyWide$PFSevent == 1 & simStudyWide$PFStimeCal <= actual2[2] &
    simStudyWide$trt == 2]), 7)
})


# censoringByNumberEvents----

test_that("censoringByNumberEvents works as expected", {
  transition1 <- weibull_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6, p01 = 0.8, p02 = 0.9, p12 = 1)
  transition2 <- weibull_transition(h01 = 1, h02 = 1.3, h12 = 1.7, p01 = 1.1, p02 = 0.9, p12 = 1.1)

  simStudy <- getOneClinicalTrial(
    nPat = c(30, 30), transitionByArm = list(transition1, transition2),
    dropout = list(rate = 0.5, time = 12),
    accrual = list(param = "time", value = 7)
  )
  simStudyWide <- getDatasetWideFormat(simStudy)

  actual <- censoringByNumberEvents(data = simStudyWide, eventNum = 12, typeEvent = "OS")
  actual2 <- censoringByNumberEvents(data = simStudyWide, eventNum = 16, typeEvent = "PFS")
  expect_equal(length(actual$id[actual$OSevent == 1]), 12)
  expect_equal(length(actual2$id[actual2$PFSevent == 1]), 16)
})

# trackEventsPerTrial----

test_that("trackEventsPerTrial works as expected by arm", {
  transition1 <- weibull_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6, p01 = 0.8, p02 = 0.9, p12 = 1)
  transition2 <- weibull_transition(h01 = 1, h02 = 1.3, h12 = 1.7, p01 = 1.1, p02 = 0.9, p12 = 1.1)

  simStudy <- getOneClinicalTrial(
    nPat = c(30, 30), transitionByArm = list(transition1, transition2),
    dropout = list(rate = 0.5, time = 2),
    accrual = list(param = "time", value = 5)
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


# getCensoredData----

test_that("getCensoredData works as expected", {
  time <- c(1, 1.2, 1.4, 1.5)
  event <- c(1, 0, 1, 1)
  data <- data.frame(id = 1:4, censTimeInd = c(1.4, 1.3, 1.5, 0.8), recruitTime = c(1, 1, 2, 3))

  actual <- getCensoredData(time = time, event = event, data = data)

  expected <- data.frame(
    time = c(1, 1.2, 1.4, 0.8), Censored = c(0, 1, 0, 1),
    event = c(1, 0, 1, 0), timeCal = c(2, 2.2, 3.4, 3.8)
  )
  expect_equal(actual, expected)
})

# getEventsAll----
test_that("getEventsAll  works as expected", {
  transition1 <- weibull_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6, p01 = 0.8, p02 = 0.9, p12 = 1)
  transition2 <- weibull_transition(h01 = 1, h02 = 1.3, h12 = 1.7, p01 = 1.1, p02 = 0.9, p12 = 1.1)
  set.seed(12345)
  simStudy <- getOneClinicalTrial(
    nPat = c(5, 5), transitionByArm = list(transition1, transition2),
    dropout = list(rate = 0.5, time = 2),
    accrual = list(param = "time", value = 5)
  )
  simStudyWide <- getDatasetWideFormat(simStudy)

  actual <- getEventsAll(data = simStudyWide, t = 3.5)
  expect_equal(actual, c(Recruited = 7, Censored = 3, UnderObs = 1))
})



# getNumberEvents----
test_that("getNumberEvents works as expected", {
  event <- c(1, 0, 1, 1, 1, 1)
  time <- c(1.2, 1.5, 7, 2.4, 3.5, 1.8)
  actual <- getNumberEvents(event = event, time = time, t = 3)
  expect_equal(actual, 3)
})
