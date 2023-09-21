# ExpHazOS ----

test_that("ExpHazOS works as expected", {
  actual <- ExpHazOS(1, 0.2, 1.1, 0.8)
  expect_equal(actual, 1.038192, tolerance = 1e-3)

  actual2 <- ExpHazOS(2, 0.4, 1.4, 0.1)
  expect_equal(actual2, 0.266345, tolerance = 1e-3)
})

test_that("ExpHazOS works also with vector of times t", {
  result <- ExpHazOS(t = c(1:3), 0.2, 1.1, 0.8)
  expected <- c(1.0381919, 0.9777975, 0.9253826)
  expect_equal(result, expected, tolerance = 1e-3)
})

# avgHRIntegExpOS ----

test_that("avgHRIntegExpOS works as expected", {
  h0 <- list(h01 = 0.18, h02 = 0.06, h12 = 0.17)
  h1 <- list(h01 = 0.23, h02 = 0.07, h12 = 0.19)
  actual <- avgHRIntegExpOS(x = 5, h01 = 0.2, h02 = 0.5, h12 = 0.7, h0 = h0, h1 = h1, alpha = 0.5)
  expect_equal(actual, 0.297362, tolerance = 1e-3)

  actual2 <- avgHRIntegExpOS(x = 1, h01 = 0.3, h02 = 0.4, h12 = 0.5, h0 = h0, h1 = h1, alpha = 0.8)
  expect_equal(actual2, 0.3764565, tolerance = 1e-3)
})

# avgHRExpOS ----

test_that("avgHRExpOS works as expected", {
  transition1 <- exponential_transition(h01 = 0.1, h02 = 0.4, h12 = 0.3)
  transition2 <- exponential_transition(h01 = 0.06, h02 = 0.3, h12 = 0.3)
  transitionList1 <- list(transition1, transition2)
  transitionList2 <- list(transition1, transition1)

  actual <- avgHRExpOS(transitionByArm = transitionList1, alpha = 0.5, upper = 100)
  expect_equal(actual, 0.8038991, tolerance = 1e-3)

  actual2 <- avgHRExpOS(transitionByArm = transitionList2, alpha = 0.5, upper = 100)
  expect_equal(actual2, 1, tolerance = 1e-3)
})
