# ExpSurvPFS ----

test_that("ExpSurvPFS works as expected", {
  actual <- ExpSurvPFS(2, 0.3, 1.8)
  expect_equal(actual, 0.01499558)

  actual2 <- ExpSurvPFS(0, 0.3, 1.8)
  expect_equal(actual2, 1)
})



# ExpSurvOS ----

test_that("ExpSurvOS works as expected", {
  actual <- ExpSurvOS(1, 0.7, 0.5, 0.8)
  expect_equal(round(actual, 5), 0.56043)

  actual2 <- ExpSurvOS(0, 0.7, 0.5, 0.8)
  expect_equal(actual2, 1)
})



# WeibSurvPFS ----
test_that("WeibSurvPFS works as expected", {
  actual <- WeibSurvPFS(2:8, 0.3, 1.8, 1, 1)
  expect_equal(actual, ExpSurvPFS(2:8, 0.3, 1.8))

  actual2 <- WeibSurvPFS(0, 0.3, 1.8, 1.2, 0.8)
  expect_equal(actual2, 1)
})


# WeibSurvOS ----
test_that("WeibSurvOS works as expected", {
  actual <- WeibSurvOS(1:5, 0.7, 0.5, 0.8, 1, 1, 1)
  expect_equal(actual, ExpSurvOS(1:5, 0.7, 0.5, 0.8))

  actual2 <- WeibSurvOS(0, 0.7, 0.5, 0.8, 1.2, 1, 0.9)
  expect_equal(actual2, 1)
})

# pwA ----
test_that("pwA works as expected", {
  actual <- pwA(c(2, 4), c(0.7, 0.9), c(0, 3))
  expect_equal(actual, c(0.7 * 2, 0.7 * 3 + 0.9))
})



# PWCsurvPFS ----
test_that("PWCsurvPFS works as expected", {
  actual <- PWCsurvPFS(1:3, c(0.7, 0.9), c(0.5, 1), c(0, 3), c(0, 7))
  expect_equal(actual, ExpSurvPFS(1:3, 0.7, 0.5))

  actual2 <- PWCsurvPFS(0, c(0.7, 0.9), c(0.5, 0.2), c(0, 9), c(0, 7))
  expect_equal(actual2, 1)

  actual3 <- PWCsurvPFS(2, c(0.7, 0.9), c(0.5, 0.2), c(0, 1), c(0, 0.8))
  expect_equal(actual3, 0.1064585)
})




# PWCsurvOS ----
test_that("PWCsurvOS works as expected", {
  actual <- PWCsurvOS(c(0, 1, 2, 1.4), c(0.7, 0.9), c(0.5, 1), c(0.9, 1.2), c(0, 4), c(0, 3), c(0, 7))
  expect_equal(actual, ExpSurvOS(c(0, 1, 2, 1.4), 0.7, 0.5, 0.9))
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
