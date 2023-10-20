# getOneToTwoRows ----

test_that("getOneToTwoRows works as expected", {
  simDataOne <- data.frame(
    id = c(1:3), to = c(1, 1, 1), from = c(0, 0, 0), entry = c(0, 0, 0),
    exit = c(0.43, 0.15, 0.01), censTime = c(0.9, 5.9, 0.5)
  )

  transition <- exponential_transition(1, 1.6, 0.3)
  actual <- getOneToTwoRows(simDataOne, transition)

  expect_true(length(actual$id) == 3)
  expect_true(all(actual$from == 1))
  expect_true(all(actual$to %in% c("cens", 2)))
  expect_true(all(actual$exit > actual$entry))
  expect_true(all(actual$entry == simDataOne$exit))
  actualCens <- actual[actual$to == "cens", ]
  expect_true(all(actualCens$censTime == actualCens$exit))
})

# getSimulatedData ----

test_that("getSimulatedData works as expected if no random censoring is present", {
  actual <- getSimulatedData(10,
    transition = exponential_transition(h01 = 1, h02 = 1.5, h12 = 1),
    dropout = list(rate = 0, time = 1),
    accrual = list(param = "time", value = 5)
  )

  # at least one row per patient and not more than two per patient.
  expect_true(length(actual$id) >= 10 & length(actual$id) <= 20)
  expect_named(actual, c("id", "from", "to", "entry", "exit", "entryAct", "exitAct", "censAct"))

  expect_true(all(actual$exit > actual$entry))
  # no censoring present?
  expect_true(all(actual$to %in% c(1, 2)))
  # illness-death model
  expect_true(all(actual$from %in% c(0, 1)))
})


test_that("getSimulatedData works as expected if random censoring is present", {
  actual <- getSimulatedData(10,
    transition = exponential_transition(h01 = 1, h02 = 1.5, h12 = 1),
    dropout = list(rate = 0.5, time = 1),
    accrual = list(param = "time", value = 5)
  )

  # at least one row per patient and not more than two per patient.
  expect_true(length(actual$id) >= 10 & length(actual$id) <= 20)
  expect_named(actual, c("id", "from", "to", "entry", "exit", "entryAct", "exitAct", "censAct"))

  expect_true(all(actual$exit > actual$entry))
  # no censoring present?
  expect_true(all(actual$to %in% c(1, 2, "cens")))
  # illness-death model
  expect_true(all(actual$from %in% c(0, 1)))
})

test_that("getSimulatedData creates expected data set for exponential transitions", {
  set.seed(1243)
  actual <- getSimulatedData(4,
    transition = exponential_transition(h01 = 1, h02 = 1.5, h12 = 1),
    dropout = list(rate = 0.5, time = 1),
    accrual = list(param = "time", value = 5)
  )
  row1 <- data.frame(
    id = 1, from = 0, to = "2", entry = 0.0000000, exit = 0.10917836889122016,
    entryAct = 3.6695390287786722, exitAct = 3.7787173976698925, censAct = 5.2614395397317368,
    stringsAsFactors = FALSE
  )
  expect_equal(actual[1, ], row1)
})


test_that("getSimulatedData creates expected data set for Weibull transitions", {
  set.seed(1243)
  actual <- getSimulatedData(4,
    transition = weibull_transition(
      h01 = 1, h02 = 1.5, h12 = 1,
      p01 = 1.1, p02 = 0.8, p12 = 1.2
    ),
    dropout = list(rate = 0.5, time = 1),
    accrual = list(param = "time", value = 5)
  )
  row1 <- data.frame(
    id = 1, from = 0, to = "2", entry = 0.0000000, exit = 0.0842096040823,
    entryAct = 3.66953902878, exitAct = 3.75374863286, censAct = 5.26143953973,
    stringsAsFactors = FALSE
  )
  expect_equal(actual[1, ], row1)
})

test_that("getSimulatedData works also without progression to death transitions", {
  set.seed(31)
  result <- expect_silent(getSimulatedData(
    N = 1L,
    transition = piecewise_exponential(
      pw01 = c(0, 2),
      pw02 = c(0, 2),
      pw12 = c(0, 2),
      h01 = rep(0.04781994, 2),
      h02 = rep(0.03882345, 2),
      h12 = rep(0.31313900, 2)
    ),
    dropout = list(rate = 0, time = 1),
    accrual = list(param = "intensity", value = 5)
  ))
  expect_data_frame(result, nrows = 1L)
  expect_false(any(result$to == 1))
  expect_false(any(result$from == 1))
})
