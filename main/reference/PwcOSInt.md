# Helper Function of `PWCsurvOS()`

Helper Function of
[`PWCsurvOS()`](https://insightsengineering.github.io/simIDM/reference/PWCsurvOS.md)

## Usage

``` r
PwcOSInt(x, t, h01, h02, h12, pw01, pw02, pw12)
```

## Arguments

- x:

  (`numeric`)\
  variable of integration.

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

Numeric results of the integrand used to calculate the OS survival
function for piecewise constant transition hazards, see `PWCsurvOS`.
