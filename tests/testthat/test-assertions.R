# assert_positive_number ----

test_that("assert_positive_number works as expected", {
  expect_silent(assert_positive_number(1.4))
  expect_silent(assert_positive_number(0, zero_ok = TRUE))
  expect_error(assert_positive_number(0))
  expect_error(assert_positive_number(-1))
  expect_error(assert_positive_number("bla"))
})
