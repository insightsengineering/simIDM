# Helper Function for `WeibSurvOS()`

Helper Function for
[`WeibSurvOS()`](https://insightsengineering.github.io/simIDM/reference/WeibSurvOS.md)

## Usage

``` r
WeibOSInteg(x, t, h01, h02, h12, p01, p02, p12)
```

## Arguments

- x:

  (`numeric`)\
  variable of integration.

- t:

  (`numeric`)\
  study time-points.

- h01:

  (positive `number`)\
  transition hazard for 0 to 1 transition.

- h02:

  (positive `number`)\
  transition hazard for 0 to 2 transition.

- h12:

  (positive `number`)\
  transition hazard for 1 to 2 transition.

- p01:

  (positive `number`)\
  rate parameter of Weibull distribution for `h01`.

- p02:

  (positive `number`)\
  rate parameter of Weibull distribution for `h02`.

- p12:

  (positive `number`)\
  rate parameter of Weibull distribution for `h12`.

## Value

Numeric results of the integrand used to calculate the OS survival
function for Weibull transition hazards, see
[`WeibSurvOS()`](https://insightsengineering.github.io/simIDM/reference/WeibSurvOS.md).
