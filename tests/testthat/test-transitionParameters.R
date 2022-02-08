# exponential_transition ----

test_that("exponential_transition works as expected", {
  actual <- exponential_transition(1, 1.2, 0.8)

  expect_identical(actual$hazards, list(h01 = 1, h02 = 1.2, h12 = 0.8))
  expect_identical(actual$family, "exponential")
  expect_identical(actual$intervals, list(pw01 = 0, pw02 = 0, pw12 = 0))
  expect_identical(actual$weibull_rates, list(p01 = 1, p02 = 1, p12 = 1))

  expect_error(exponential_transition(-1, 1.2, 0.8))
  expect_error(exponential_transition(1, 1.2, 0.8, 6))
})

# piecewise_exponential ----

test_that("piecewise_exponential works as expected", {
  actual <- piecewise_exponential(
    h01 = c(1, 1, 1), h02 = c(1.5, 0.5, 1), h12 = c(1, 1, 1),
    pw01 = c(0, 3, 8), pw02 = c(0, 6, 7), pw12 = c(0, 8, 9)
  )

  expect_identical(actual$hazards, list(h01 = c(1, 1, 1), h02 = c(1.5, 0.5, 1), h12 = c(1, 1, 1)))
  expect_identical(actual$family, "piecewise exponential")
  expect_identical(actual$intervals, list(pw01 = c(0, 3, 8), pw02 = c(0, 6, 7), pw12 = c(0, 8, 9)))
  expect_identical(actual$weibull_rates, list(p01 = 1, p02 = 1, p12 = 1))

  # first element piecewise has to be 0
  expect_error(piecewise_exponential(
    h01 = c(1, 1, 1), h02 = c(1.5, 0.5, 1), h12 = c(1, 1, 1),
    pw01 = c(4, 3, 8), pw02 = c(0, 6, 7), pw12 = c(0, 8, 9)
  ))
  ## hazards must be positive numbers
  expect_error(piecewise_exponential(
    h01 = c(-1, 1, 1), h02 = c(1.5, 0.5, 1), h12 = c(1, 1, 1),
    pw01 = c(0, 3, 8), pw02 = c(0, 6, 7), pw12 = c(0, 8, 9)
  ))
  # intervals must be increasing
  expect_error(piecewise_exponential(
    h01 = c(1, 1, 1), h02 = c(1.5, 0.5, 1), h12 = c(1, 1, 1),
    pw01 = c(0, 3, 8), pw02 = c(0, 6, 3), pw12 = c(0, 8, 9)
  ))
})

# weibull_transition ----

test_that("weibull_transition works as expected", {
  actual <- weibull_transition(
    h01 = 1, h02 = 1.3, h12 = 0.5,
    p01 = 1.2, p02 = 1.3, p12 = 0.5
  )


  expect_identical(actual$hazards, list(h01 = 1, h02 = 1.3, h12 = 0.5))
  expect_identical(actual$family, "Weibull")
  expect_identical(actual$intervals, list(pw01 = 0, pw02 = 0, pw12 = 0))
  expect_identical(actual$weibull_rates, list(p01 = 1.2, p02 = 1.3, p12 = 0.5))

  expect_error(weibull_transition(
    h01 = 1, h02 = -1.3, h12 = 0.5,
    p01 = 1.2, p02 = 1.3, p12 = 0.5
  ))
  expect_error(weibull_transition(
    h01 = 1, h02 = 1.3, h12 = 0.5,
    p01 = 1.2, p02 = 1.3, p12 = -1
  ))
})
