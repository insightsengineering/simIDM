# Assertion for vector describing intervals

We define an intervals vector to always start with 0, and contain unique
ordered time points.

## Usage

``` r
assert_intervals(x, y)
```

## Arguments

- x:

  what to check.

- y:

  (`count`)\
  required length of `y`.

## Value

Raises an error if `x` is not an intervals vector starting with 0.

## Examples

``` r
assert_intervals(c(0, 5, 7), 3)
```
