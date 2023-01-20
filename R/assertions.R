#' Assertion for Positive Number
#'
#' @param x what to check.
#' @param zero_ok (`flag`)\cr whether `x` can be zero or not.
#'
#' @return Raises an error if `x` is not a single positive (or non-negative) number.
#' @export
#'
#' @examples
#' assert_positive_number(3.2)
#' assert_positive_number(0, zero_ok = TRUE)
assert_positive_number <- function(x,
                                   zero_ok = FALSE) {
  assert_number(x)
  assert_flag(zero_ok)
  if (zero_ok) {
    assert_true(x >= 0)
  } else {
    assert_true(x > 0)
  }
}

#' Assertion for vector describing intervals
#'
#' We define an intervals vector to always start with 0, and contain
#' unique ordered time points.
#'
#' @param x what to check.
#' @param y (`count`)\cr required length of `y`.
#'
#' @return Raises an error if `x` is not an intervals vector starting with 0.
#' @export
#'
#' @examples
#' assert_intervals(c(0, 5, 7), 3)
assert_intervals <- function(x, y) {
  assert_numeric(
    x,
    lower = 0,
    any.missing = FALSE,
    all.missing = FALSE,
    len = y,
    unique = TRUE,
    sorted = TRUE
  )
  assert_true(x[1] == 0)
}
