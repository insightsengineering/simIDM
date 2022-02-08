# assert_positive_number ----

test_that("assert_positive_number works as expected", {
  expect_silent(assert_positive_number(1.4))
  expect_silent(assert_positive_number(0, zero_ok = TRUE))
  expect_error(assert_positive_number(0))
  expect_error(assert_positive_number(-1))
  expect_error(assert_positive_number("bla"))
})

# assert_intervals ----

test_that("assert_intervals works as expected", {
  expect_silent(assert_intervals(c(0, 1, 3, 4), 4))
  expect_silent(assert_intervals(c(0, 5.5, 7.7), 3))
  expect_error(assert_intervals(c(0, 5.5, 7.7), 2))
  expect_error(assert_intervals(c(1, 5.5, 7.7), 3))
  expect_error(assert_intervals(c(0, 7, 7), 3))
  expect_error(assert_intervals(c(0, 7, 2), 3))
})
