# getPCWHazard ----

test_that("getPCWHazard works as expected", {
  actual <- getPCWHazard(c(0.8, 1.1, 1), c(0, 5, 8), c(1, 5, 7, 9))
  expect_identical(actual, c(0.8, 1.1, 1.1, 1))
})
