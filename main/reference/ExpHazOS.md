# OS Hazard Function from Constant Transition Hazards

OS Hazard Function from Constant Transition Hazards

## Usage

``` r
ExpHazOS(t, h01, h02, h12)
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

This returns the value of the OS hazard function at time t.

## Examples

``` r
ExpHazOS(c(1:5), 0.2, 1.1, 0.8)
#> [1] 1.0381919 0.9777975 0.9253826 0.8843734 0.8548146
```
