# OS Survival Function from Constant Transition Hazards

OS Survival Function from Constant Transition Hazards

## Usage

``` r
ExpSurvOS(t, h01, h02, h12)
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

## Value

This returns the value of OS survival function at time t.

## Examples

``` r
ExpSurvOS(c(1:5), 0.2, 0.4, 0.1)
#> [1] 0.6912219 0.5082088 0.3955066 0.3225588 0.2724845
```
