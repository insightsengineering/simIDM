# ExpSurvPFS ----

test_that("ExpSurvPFS works as expected", {
  actual <- ExpSurvPFS(2, 0.3, 1.8)
  expect_equal(actual, 0.01499558, tolerance = 1e-3)

  actual2 <- ExpSurvPFS(0, 0.3, 1.8)
  expect_equal(actual2, 1)
})

test_that("ExpSurvPFS works also with vector of times t", {
  result <- ExpSurvPFS(t = c(1, 0.5, 3), 0.3, 1.8)
  expected <- c(0.1225, 0.3499, 0.0018)
  expect_equal(result, expected, tolerance = 1e-3)
})

# ExpSurvOS ----

test_that("ExpSurvOS works as expected", {
  actual <- ExpSurvOS(1, 0.7, 0.5, 0.8)
  expect_equal(round(actual, 5), 0.56043)

  actual2 <- ExpSurvOS(0, 0.7, 0.5, 0.8)
  expect_equal(actual2, 1)

  actual3 <- ExpSurvOS(1000, 1.2, 1.5, 1.6)
  expect_equal(actual3, 0)
})

test_that("ExpSurvOS works also with vector of times t", {
  result <- ExpSurvOS(t = c(1, 0.5, 3), 0.7, 0.5, 0.8)
  expected <- c(0.5604, 0.7615, 0.1383)
  expect_equal(result, expected, tolerance = 1e-3)
})

# WeibSurvPFS ----

test_that("WeibSurvPFS works as expected", {
  actual <- WeibSurvPFS(2:8, 0.3, 1.8, 1, 1)
  expect_equal(actual, ExpSurvPFS(2:8, 0.3, 1.8))

  actual2 <- WeibSurvPFS(0, 0.3, 1.8, 1.2, 0.8)
  expect_equal(actual2, 1)
})

# integrateVector ----

test_that("integrateVector works as expected", {
  integrand <- function(x) x^2
  upper <- c(1, 0.4, 1)

  actual <- integrateVector(integrand, upper = upper)
  expected <- c(
    integrate(integrand, 0, 1)$value,
    integrate(integrand, 0, 0.4)$value,
    integrate(integrand, 0, 1)$value
  )
  expect_equal(actual, expected)
})

# WeibOSInteg ----

test_that("WeibOSInteg works as expected with scalar x", {
  result <- expect_silent(WeibOSInteg(4, 2:6, 0.2, 0.5, 2.1, 1.2, 0.9, 1))
  expected <- c(5.368515, 0.657409, 0.080504, 0.009858, 0.001207)
  expect_equal(result, expected, tolerance = 1e-4)
})

test_that("WeibOSInteg works as expected with vector x", {
  result <- expect_silent(WeibOSInteg(1:5, 2:6, 0.2, 0.5, 2.1, 1.2, 0.9, 1))
  expected <- c(0.06081, 0.034948, 0.018842, 0.009858, 0.005061)
  expect_equal(result, expected, tolerance = 1e-4)
})

test_that("WeibOSInteg works as expected with scalar t", {
  result <- expect_silent(WeibOSInteg(2:6, 4, 0.2, 0.5, 2.1, 1.2, 0.9, 1))
  expected <- c(0.00428, 0.018842, 0.080504, 0.337499, 1.395583)
  expect_equal(result, expected, tolerance = 1e-4)
})

# WeibSurvOS ----

test_that("WeibSurvOS works as expected", {
  actual <- WeibSurvOS(1:5, 0.7, 0.5, 0.8, 1, 1, 1)
  expect_equal(actual, ExpSurvOS(1:5, 0.7, 0.5, 0.8))

  actual2 <- WeibSurvOS(0, 0.7, 0.5, 0.8, 1.2, 1, 0.9)
  expect_equal(actual2, 1)

  actual3 <- WeibSurvOS(1000, 1, 1.1, 1.2, 1.3, 0.8, 1.4)
  expect_equal(actual3, 0)
})

# pwA ----

test_that("pwA works as expected", {
  actual <- pwA(c(2, 4), c(0.7, 0.9), c(0, 3))
  expect_equal(actual, c(0.7 * 2, 0.7 * 3 + 0.9))
})

# PWCsurvPFS ----

test_that("PWCsurvPFS works as expected", {
  actual <- PWCsurvPFS(1:3, c(0.7, 0.9), c(0.5, 1), c(0, 3), c(0, 7))
  expect_equal(actual, ExpSurvPFS(1:3, 0.7, 0.5), tolerance = 1e-3)

  actual2 <- PWCsurvPFS(0, c(0.7, 0.9), c(0.5, 0.2), c(0, 9), c(0, 7))
  expect_equal(actual2, 1)

  actual3 <- PWCsurvPFS(2, c(0.7, 0.9), c(0.5, 0.2), c(0, 1), c(0, 0.8))
  expect_equal(actual3, 0.1064585, tolerance = 1e-3)
})

# PWCsurvOS ----

test_that("PWCsurvOS works as expected", {
  actual <- PWCsurvOS(c(0, 1, 2, 1.4), c(0.7, 0.9), c(0.5, 1), c(0.9, 1.2), c(0, 4), c(0, 3), c(0, 7))
  expect_equal(actual, ExpSurvOS(c(0, 1, 2, 1.4), 0.7, 0.5, 0.9))

  actual2 <- PWCsurvOS(100, 0.06, 1, 8, 0, 0, 0)
  expect_equal(actual2, 0, tolerance = 1e-10)

  actual3 <- PWCsurvOS(2000, c(0.7, 0.9), c(0.5, 1), c(0.9, 14), c(0, 4), c(0, 3), c(0, 7))
  expect_equal(actual3, 0, tolerance = 1e-10)
})

test_that("PWCsurvOS does not return values larger than 1", {
  result <- PWCsurvOS(3.13, c(0, 1, 1), c(0, 0.5, 1), c(0, 1, 1), c(0, 3, 8), c(0, 6, 7), c(0, 8, 9))
  expect_lte(result, 1)
})

test_that("PWCsurvOS is monotonically decreasing and not larger than 1", {
  b1 <- 3
  b2 <- 6
  t_values <- seq(2.9, 6.2, length.out = 100)
  assert_true(all(diff(t_values) > 0))
  result <- PWCsurvOS(t_values, c(0, 1, 1), c(0, 0.5, 1), c(0, 1, 1), c(0, b1, 8), c(0, b2, 7), c(0, 8, 9))
  expect_true(all(diff(result) <= 0))
  expect_true(all(result <= 1))
  expect_true(all(result >= 0))
})

test_that("PWCsurvOS gives equal results as the numerical integration", {
  t <- c(0.1, 1, 5.2, 3, 10.5, 2, 7.3, 1.5, 5.2, 60.5, 0.3, 20.6, 0.01)
  h01 <- c(1, 0.2, 3.2)
  pw01 <- c(0, 2.5, 15.3)
  h02 <- 4
  pw02 <- 0
  h12 <- c(5.1, 4.2, 3.1, 7.3, 50.5)
  pw12 <- c(0, 7.3, 20, 50.2, 70.5)

  result <- PWCsurvOS(t, h01, h02, h12, pw01, pw02, pw12)
  expected <- PWCsurvPFS(t, h01, h02, pw01, pw02) +
    sapply(t, function(t) {
      integrateVector(PwcOSInt,
        upper = t,
        t = t,
        h01 = h01,
        h02 = h02,
        h12 = h12,
        pw01 = pw01,
        pw02 = pw02,
        pw12 = pw12
      )
    })
  expect_equal(result, expected, tolerance = 1e-6)
})

# PwcOSInt ----

test_that("PwcOSInt works as expected", {
  actual <- PwcOSInt(1, 1.1, c(0.3, 0.5), c(0.5, 0.8), c(0.7, 1), c(0, 4), c(0, 8), c(0, 3))
  expect_equal(actual, 0.12568546, tolerance = 1e-6)

  actual2 <- PwcOSInt(4, 4, c(0.3, 0.5), c(0.5, 0.8), c(0.7, 100), c(0, 2), c(0, 3), c(0, 3))
  expect_equal(actual2, 0.010120956, tolerance = 1e-6)

  actual3 <- PwcOSInt(700, 700, c(0.3, 0.5), c(0.5, 0.8), c(0.7, 100), c(0, 2), c(0, 3), c(0, 3))
  expect_equal(actual3, 0, tolerance = 1e-6)
})

# singleExpQuantOS ----

test_that("singleExpQuantOS works as expected", {
  actual <- singleExpQuantOS(1 / 2, 0.2, 0.5, 2.1)
  expect_equal(actual, 1.144539, tolerance = 1e-5)
})

# ExpQuantOS ----

test_that("ExpQuantOS works as expected", {
  actual <- ExpQuantOS(1 / 2, 0.2, 0.5, 2.1)
  expect_equal(actual, 1.144539, tolerance = 1e-5)
})

test_that("ExpQuantOS works also with vector of quantiles q", {
  actual <- ExpQuantOS(q = c(0.1, 0.3, 0.7), 0.2, 0.5, 2.1)
  expected <- c(3.4788, 1.8981, 0.6237)
  expect_equal(actual, expected, tolerance = 1e-5)
})
