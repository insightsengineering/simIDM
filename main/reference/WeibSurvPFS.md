# PFS Survival Function from Weibull Transition Hazards

PFS Survival Function from Weibull Transition Hazards

## Usage

``` r
WeibSurvPFS(t, h01, h02, p01, p02)
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

- p01:

  (positive `number`)\
  rate parameter of Weibull distribution for `h01`.

- p02:

  (positive `number`)\
  rate parameter of Weibull distribution for `h02`.

## Value

This returns the value of PFS survival function at time t.

## Examples

``` r
WeibSurvPFS(c(1:5), 0.2, 0.5, 1.2, 0.9)
#> [1] 0.49658530 0.24845033 0.12351703 0.06101061 0.02995439
```
