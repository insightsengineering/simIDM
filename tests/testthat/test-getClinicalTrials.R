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
    id = 1, from = 0, to = "1", entry = 0.0000000, exit = 0.1953989333804476,
    entryAct = 0, exitAct = 0.1953989333804476, censAct = 19.10280613143678, trt = 1,
    stringsAsFactors = FALSE
  )
  expect_equal(actual[1, ], row1)
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
