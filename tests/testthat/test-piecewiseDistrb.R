# getPCWdistr ----
test_that("getPCWDistr works as expected", {
  U <- runif(1)
  actual <- getPCWDistr(U, c(1.1, 0.8), c(0, 5), 6)
  expected <- -log(1 - U) / 0.8
  expect_equal(actual, expected)

  set.seed(1225)
  U2 <- runif(2)
  actual2 <- getPCWDistr(U2, c(1.1, 0.5, 0.4), c(0, 0.6, 5), c(0, 4.2))

  # 1-survival function of the actual time-points.
  expected2 <- 1 - exp(-(1.1 * 0.6) - (actual2[1] - 0.6) * 0.5)

  expected3 <- 1 - exp(-((actual2[2] + 4.2) - 4.2) * 0.5)

  expect_equal(U2[1], expected2)
  expect_equal(U2[2], expected3)
})
