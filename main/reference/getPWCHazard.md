# Piecewise Constant Hazard Values

This returns piecewise constant hazard values at specified time points.

## Usage

``` r
getPWCHazard(haz, pw, x)
```

## Arguments

- haz:

  (`numeric`)\
  piecewise constant input hazard.

- pw:

  (`numeric`)\
  time intervals for the piecewise constant hazard.

- x:

  (`numeric`)\
  time-points.

## Value

Hazard values at input time-points.

## Examples

``` r
getPWCHazard(c(1, 1.2, 1.4), c(0, 2, 3), c(1, 4, 6))
#> [1] 1.0 1.4 1.4
```
