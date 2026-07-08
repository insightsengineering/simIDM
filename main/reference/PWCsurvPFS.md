# PFS Survival Function from Piecewise Constant Hazards

PFS Survival Function from Piecewise Constant Hazards

## Usage

``` r
PWCsurvPFS(t, h01, h02, pw01, pw02)
```

## Arguments

- t:

  (`numeric`)\
  study time-points.

- h01:

  (`numeric vector`)\
  constant transition hazards for 0 to 1 transition.

- h02:

  (`numeric vector`)\
  constant transition hazards for 0 to 2 transition.

- pw01:

  (`numeric vector`)\
  time intervals for the piecewise constant hazards `h01`.

- pw02:

  (`numeric vector`)\
  time intervals for the piecewise constant hazards `h02`.

## Value

This returns the value of PFS survival function at time t.

## Examples

``` r
PWCsurvPFS(1:5, c(0.3, 0.5), c(0.5, 0.8), c(0, 4), c(0, 8))
#> [1] 0.44932896 0.20189652 0.09071795 0.04076220 0.01499558
```
