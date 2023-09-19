# ExpHazOS ----

test_that("ExpHazOS works as expected", {
  actual <- ExpHazOS(1, 0.2, 1.1, 0.8)
  expect_equal(actual, 1.038192, tolerance = 1e-3)

  actual2 <- ExpHazOS(2, 0.4, 1.4, 0)
  expect_equal(actual2, 0.1571141, tolerance = 1e-3)
})

test_that("ExpHazOS works also with vector of times t", {
  result <- ExpHazOS(t = c(1:3), 0.2, 1.1, 0.8)
  expected <- c(1.0381919, 0.9777975, 0.9253826)
  expect_equal(result, expected, tolerance = 1e-3)
})
