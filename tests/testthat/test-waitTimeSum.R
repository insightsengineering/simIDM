# getWaitTimeSum ----
test_that("getWaitTimeSum  works as expected for sum of constant hazards", {
  U <- 0.65
  actual <- getWaitTimeSum(U, haz1 = 0.8, haz2 = 1.2, p1 = 1, p2 = 1, entry = 0)
  expected <- -log(1 - U) / (0.8 + 1.2)
  expect_equal(actual, expected)
})


test_that("getWaitTimeSum  works as expected for sum of hazards (entry = 0)", {
  U <- c(0.23, 0.86)
  actual <- getWaitTimeSum(U, haz1 = 0.8, haz2 = 1.2, p1 = 1.5, p2 = 0.7, entry = c(0, 0))
  expected <- c((0.8 * actual[1])^1.5 + (1.2 * actual[1])^0.7, (0.8 * actual[2])^1.5 + (1.2 * actual[2])^0.7)
  expect_equal(-log(1 - U), expected)
})

test_that("getWaitTimeSum  works as expected for sum of hazards (entry != 0)", {
  U <- 0.4
  actual <- getWaitTimeSum(U, haz1 = 0.8, haz2 = 1.2, p1 = 1.5, p2 = 0.7, entry = 5)
  expected <- (0.8 * (actual + 5))^1.5 - (0.8 * 5)^1.5 - (1.2 * 5)^0.7 + (1.2 * (actual + 5))^0.7
  expect_equal(-log(1 - U), expected)
})
