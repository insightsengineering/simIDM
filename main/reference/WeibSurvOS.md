# OS Survival Function from Weibull Transition Hazards

OS Survival Function from Weibull Transition Hazards

## Usage

``` r
WeibSurvOS(t, h01, h02, h12, p01, p02, p12)
```

## Arguments

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

This returns the value of OS survival function at time t.

## Examples

``` r
WeibSurvOS(c(1:5), 0.2, 0.5, 2.1, 1.2, 0.9, 1)
#> [1] 0.55296799 0.29049137 0.14798085 0.07421549 0.03684786
```
