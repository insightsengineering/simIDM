# Event Times Distributed as Sum of Weibull

This returns event times with a distribution resulting from the sum of
two Weibull distributed random variables using the inversion method.

## Usage

``` r
getWaitTimeSum(U, haz1, haz2, p1, p2, entry)
```

## Arguments

- U:

  (`numeric`)\
  uniformly distributed random variables.

- haz1:

  (positive `number`)\
  first summand (constant hazard).

- haz2:

  (positive `number`)\
  second summand (constant hazard).

- p1:

  (positive `number`)\
  rate parameter of Weibull distribution for `haz1`.

- p2:

  (positive `number`)\
  rate parameter of Weibull distribution for `haz2`.

- entry:

  (`numeric`)\
  the entry times in the current state.

## Value

This returns a vector with event times.

## Examples

``` r
getWaitTimeSum(U = c(0.4, 0.5), haz1 = 0.8, haz2 = 1, p1 = 1.1, p2 = 1.5, entry = c(0, 0))
#> [1] 0.3803563 0.4820279
```
