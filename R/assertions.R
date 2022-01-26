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
