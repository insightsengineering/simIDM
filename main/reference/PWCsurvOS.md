# OS Survival Function from Piecewise Constant Hazards

OS Survival Function from Piecewise Constant Hazards

## Usage

``` r
PWCsurvOS(t, h01, h02, h12, pw01, pw02, pw12)
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

- h12:

  (`numeric vector`)\
  constant transition hazards for 1 to 2 transition.

- pw01:

  (`numeric vector`)\
  time intervals for the piecewise constant hazards `h01`.

- pw02:

  (`numeric vector`)\
  time intervals for the piecewise constant hazards `h02`.

- pw12:

  (`numeric vector`)\
  time intervals for the piecewise constant hazards `h12`.

## Value

This returns the value of OS survival function at time t.

## Examples

``` r
PWCsurvOS(1:5, c(0.3, 0.5), c(0.5, 0.8), c(0.7, 1), c(0, 4), c(0, 8), c(0, 3))
#> [1] 0.59109798 0.33599786 0.18593338 0.08687340 0.03945673
```
