# Assertion for Positive Number

Assertion for Positive Number

## Usage

``` r
assert_positive_number(x, zero_ok = FALSE)
```

## Arguments

- x:

  what to check.

- zero_ok:

  (`flag`)\
  whether `x` can be zero or not.

## Value

Raises an error if `x` is not a single positive (or non-negative)
number.

## Examples

``` r
assert_positive_number(3.2)
assert_positive_number(0, zero_ok = TRUE)
```
