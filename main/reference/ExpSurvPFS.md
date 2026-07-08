# PFS Survival Function from Constant Transition Hazards

PFS Survival Function from Constant Transition Hazards

## Usage

``` r
ExpSurvPFS(t, h01, h02)
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

## Value

This returns the value of PFS survival function at time t.

## Examples

``` r
ExpSurvPFS(c(1:5), 0.2, 0.4)
#> [1] 0.54881164 0.30119421 0.16529889 0.09071795 0.04978707
```
