# survPFS ----

test_that("survPFS works as expected for exponential transition hazards", {
  transition <- exponential_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6)
  t <- c(0.4, 1.4)
  result <- expect_silent(survPFS(transition, t))
  expected <- ExpSurvPFS(t = t, h01 = transition$hazards$h01, h02 = transition$hazards$h02)
  expect_identical(result, expected)
})

test_that("survPFS works as expected for Weibull transition hazards", {
  transition <- weibull_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6, p01 = 2, p02 = 2.5, p12 = 3)
  t <- c(0.4, 1.4)
  result <- expect_silent(survPFS(transition, t))
  expected <- WeibSurvPFS(
    t = t, h01 = transition$hazards$h01, h02 = transition$hazards$h02,
    p01 = transition$weibull_rates$p01, p02 = transition$weibull_rates$p02
  )
  expect_identical(result, expected)
})

test_that("survPFS works as expected for piecewise exponential transition hazards", {
  transition <- piecewise_exponential(
    h01 = c(1, 1, 1), h02 = c(1.5, 0.5, 1), h12 = c(1, 1, 1),
    pw01 = c(0, 3, 8), pw02 = c(0, 6, 7), pw12 = c(0, 8, 9)
  )
  t <- c(0.4, 1.4)
  result <- expect_silent(survPFS(transition, t))
  expected <- PWCsurvPFS(
    t = t, h01 = transition$hazards$h01, h02 = transition$hazards$h02,
    pw01 = transition$intervals$pw01, pw02 = transition$intervals$pw02
  )
  expect_identical(result, expected)
})

# survOS ----

test_that("survOS works as expected for exponential transition hazards", {
  transition <- exponential_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6)
  t <- c(0.4, 1.4)
  result <- expect_silent(survOS(transition, t))
  expected <- ExpSurvOS(
    t = t, h01 = transition$hazards$h01, h02 = transition$hazards$h02,
    h12 = transition$hazards$h12
  )
  expect_identical(result, expected)
})

test_that("survOS works as expected for Weibull transition hazards", {
  transition <- weibull_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6, p01 = 2, p02 = 2.5, p12 = 3)
  t <- c(0.4, 1.4)
  result <- expect_silent(survOS(transition, t))
  expected <- WeibSurvOS(
    t = t, h01 = transition$hazards$h01, h02 = transition$hazards$h02,
    h12 = transition$hazards$h12, p01 = transition$weibull_rates$p01,
    p02 = transition$weibull_rates$p02, p12 = transition$weibull_rates$p12
  )
  expect_identical(result, expected)
})

test_that("survOS works as expected for piecewise exponential transition hazards", {
  transition <- piecewise_exponential(
    h01 = c(1, 1, 1), h02 = c(1.5, 0.5, 1), h12 = c(1, 1, 1),
    pw01 = c(0, 3, 8), pw02 = c(0, 6, 7), pw12 = c(0, 8, 9)
  )
  t <- c(0.4, 1.4)
  result <- expect_silent(survOS(transition, t))
  expected <- PWCsurvOS(
    t = t, h01 = transition$hazards$h01, h02 = transition$hazards$h02,
    h12 = transition$hazards$h12, pw01 = transition$intervals$pw01,
    pw02 = transition$intervals$pw02, pw12 = transition$intervals$pw12
  )
  expect_identical(result, expected)
})

# expval ----

# expvalPFSInteg ----

test_that("expvalPFSInteg works as expected", {
  transition <- exponential_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6)
  x <- c(0.4, 10.3)
  result <- expect_silent(expvalPFSInteg(x, transition))
  expected <- x * survPFS(transition, x)
  expect_identical(result, expected)
})

# expvalOSInteg ----

test_that("expvalOSInteg works as expected", {
  transition <- exponential_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6)
  x <- c(0.4, 10.3)
  result <- expect_silent(expvalOSInteg(x, transition))
  expected <- x * survOS(transition, x)
  expect_identical(result, expected)
})

# p11 ----

test_that("p11Integ works as expected", {
  transition <- exponential_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6)
  x <- c(1.45, 0.6)
  result <- expect_silent(p11Integ(x, transition))
  expected <- haz(transition = transition, t = x, trans = 3)
  expect_identical(result, expected)
})

test_that("log_p11 works as expected", {
  transition <- exponential_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6)
  s <- 1
  t <- 3
  result <- expect_silent(log_p11(transition, s, t))
  expected <- -transition$hazards$h12 * (t - s)
  expect_equal(result, expected)
})

test_that("log_p11 works as expected for multiple time points", {
  transition <- exponential_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6)
  s <- c(1, 2.2)
  t <- c(4.2, 3.2)
  result <- expect_silent(log_p11(transition, s, t))
  expected <- -transition$hazards$h12 * (t - s)
  expect_equal(result, expected)
})

# PFSOS ----

test_that("PFSOSInteg works as expected", {
  transition <- exponential_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6)
  result <- expect_silent(PFSOSInteg(1, 2, transition))
  expected <- 0.0162823
  expect_equal(result, expected, tolerance = 1e-4)
})

test_that("PFSOSInteg works as expected for vectors", {
  transition <- exponential_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6)
  result <- expect_silent(PFSOSInteg(c(1, 2), c(5, 7), transition))
  expected <- c(0.000134, 0.000492)
  expect_equal(result, expected, tolerance = 1e-3)
})

test_that("survPFSOS works as expected for exponential transition", {
  transition <- exponential_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6)
  result <- expect_silent(survPFSOS(c(0.2, 0.1, 3.5), transition))
  expected <- c(0.38618, 0.52293, 0.01085)
  expect_equal(result, expected, tolerance = 1e-4)
})

test_that("survPFSOS works as expected for Weibull transitions", {
  transition <- weibull_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6, p01 = 2, p02 = 2.5, p12 = 3)
  result <- expect_silent(survPFSOS(c(0.2, 0.1, 2.5), transition))
  expected <- c(0.75197, 0.89561, 0.00056)
  expect_equal(result, expected, tolerance = 1e-4)
})

test_that("survPFSOS works as expected for piecewise constant hazard transitions", {
  transition <- piecewise_exponential(
    h01 = c(1, 1, 1), h02 = c(1.5, 0.5, 1), h12 = c(1, 1, 1),
    pw01 = c(0, 3, 8), pw02 = c(0, 6, 7), pw12 = c(0, 8, 9)
  )
  result <- expect_silent(survPFSOS(c(0.2, 0.1, 2.5), transition))
  expected <- c(0.43021, 0.5602, 0.03678)
  expect_equal(result, expected, tolerance = 1e-4)
})

# correlation ----

test_that("corTrans works as expected", {
  transition <- weibull_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6, p01 = 2, p02 = 2.5, p12 = 3)
  result <- expect_silent(corTrans(transition))
  expected <- 0.5993243
  expect_equal(result, expected, tolerance = 1e-4)
})

test_that("corPFSOS returns correct central estimate", {
  transition <- exponential_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6)
  data <- getClinicalTrials(
    nRep = 1, nPat = c(20000), seed = 1234, datType = "1rowTransition",
    transitionByArm = list(transition), dropout = list(rate = 0.5, time = 12),
    accrual = list(param = "intensity", value = 7)
  )[[1]]
  actual <- corPFSOS(data, transition = exponential_transition(), bootstrap = FALSE)
  expect_equal(actual[[1]], 0.5764501, tolerance = 1e-3)
})
