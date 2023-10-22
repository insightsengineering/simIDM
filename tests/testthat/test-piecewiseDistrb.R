# getPCWdistr ----

test_that("getPCWDistr works as expected for constant hazards after time shift", {
  U <- 0.65
  actual <- getPCWDistr(U, c(1.1, 0.8), c(0, 5), 6)
  expected <- -log(1 - U) / 0.8
  expect_equal(actual, expected)
})

test_that("getPCWDistr works as expected for pwc hazards without time shift", {
  U <- 0.91
  actual <- getPCWDistr(U = U, haz = c(1.1, 0.5, 0.4), pw = c(0, 0.6, 5), t_0 = 0)
  # 1-survival function of the actual time-points.
  expected <- 1 - exp(-(1.1 * 0.6) - (actual - 0.6) * 0.5)
  expect_equal(U, expected)
})


test_that("getPCWDistr works as expected for piecewise constant hazards with time shift", {
  U <- c(0.45, 0.8)
  actual <- getPCWDistr(U = U, haz = c(1.1, 0.5, 0.4), pw = c(0, 0.6, 5), t_0 = c(4.2, 4.2))
  # 1-survival function of the actual time-points.
  expected <- 1 - exp(-(0.8 * 0.5) - (actual + 4.2 - 5) * 0.4)
  expect_equal(U, expected)
})

# PCWInversionMethod ----

test_that("PCWInversionMethod works as expected", {
  U <- 0.4
  logU <- log(1 - U)
  actual <- PCWInversionMethod(haz = c(1.1, 0.5, 0.4), pw = c(0, 0.6, 5), LogU = logU)
  # 1-survival function of the actual time-point.
  expected <- 1 - exp(-(actual * 1.1))
  expect_equal(U, expected)
})
