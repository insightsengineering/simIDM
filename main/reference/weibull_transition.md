# Transition Hazards for Weibull Distributed Event Times

This creates a list with class `TransitionParameters` containing
hazards, time intervals and Weibull rates for Weibull distributed event
times in an illness-death model.

## Usage

``` r
weibull_transition(h01 = 1, h02 = 1, h12 = 1, p01 = 1, p02 = 1, p12 = 1)
```

## Arguments

- h01:

  (positive `number`)\
  transition hazard for 0 to 1 transition

- h02:

  (positive `number`)\
  transition hazard for 0 to 2 transition

- h12:

  (positive `number`)\
  transition hazard for 1 to 2 transition

- p01:

  (positive `number`)\
  rate parameter of Weibull distribution for `h01`

- p02:

  (positive `number`)\
  rate parameter of Weibull distribution for `h02`

- p12:

  (positive `number`)\
  rate parameter of Weibull distribution for `h12`

## Value

List with elements `hazards`, `intervals`, `weibull_rates` and `family`
(Weibull).

## Examples

``` r
weibull_transition(h01 = 1, h02 = 1.3, h12 = 0.5, p01 = 1.2, p02 = 1.3, p12 = 0.5)
#> $hazards
#> $hazards$h01
#> [1] 1
#> 
#> $hazards$h02
#> [1] 1.3
#> 
#> $hazards$h12
#> [1] 0.5
#> 
#> 
#> $weibull_rates
#> $weibull_rates$p01
#> [1] 1.2
#> 
#> $weibull_rates$p02
#> [1] 1.3
#> 
#> $weibull_rates$p12
#> [1] 0.5
#> 
#> 
#> $intervals
#> $intervals$pw01
#> [1] 0
#> 
#> $intervals$pw02
#> [1] 0
#> 
#> $intervals$pw12
#> [1] 0
#> 
#> 
#> $family
#> [1] "Weibull"
#> 
#> attr(,"class")
#> [1] "WeibullTransition"    "TransitionParameters"
```
