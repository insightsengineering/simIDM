# Cumulative Hazard for Piecewise Constant Hazards

Cumulative Hazard for Piecewise Constant Hazards

## Usage

``` r
pwA(t, haz, pw)
```

## Arguments

- t:

  (`numeric`)\
  study time-points.

- haz:

  (`numeric vector`)\
  constant transition hazards.

- pw:

  (`numeric vector`)\
  time intervals for the piecewise constant hazards.

## Value

This returns the value of cumulative hazard at time t.

## Examples

``` r
pwA(1:5, c(0.5, 0.9), c(0, 4))
#> [1] 0.5 1.0 1.5 2.0 2.9
```
