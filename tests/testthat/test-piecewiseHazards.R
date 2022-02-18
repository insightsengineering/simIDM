# getPCWHazard ----

test_that("getPCWHazard works as expected", {
  actual <- getPCWHazard(c(0.8, 1.1, 1), c(0, 5, 8), c(1, 5, 7, 9))
  expect_identical(actual, c(0.8, 1.1, 1.1, 1))
})


# getSumPCW ----
test_that("getSumPCW works as expected", {
  actual <- getSumPCW(c(0.8, 1.1, 1), c(0.8, 1.1, 1, 0.4), c(0, 5, 8), c(0, 3, 7, 9))
  expect_equal(actual, list(hazards = c(1.6, 1.9, 2.2, 2.1, 2, 1.4), intervals = c(0, 3, 5, 7, 8, 9)))
})
