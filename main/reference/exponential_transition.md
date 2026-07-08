# Transition Hazards for Exponential Event Times

This creates a list with class `TransitionParameters` containing
hazards, time intervals and Weibull rates for exponential event times in
an illness-death model.

## Usage

``` r
exponential_transition(h01 = 1, h02 = 1, h12 = 1)
```

## Arguments

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

List with elements `hazards`, `intervals`, `weibull_rates` and `family`
(exponential).

## Examples

``` r
exponential_transition(1, 1.6, 0.3)
#> $hazards
#> $hazards$h01
#> [1] 1
#> 
#> $hazards$h02
#> [1] 1.6
#> 
#> $hazards$h12
#> [1] 0.3
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
#> $weibull_rates
#> $weibull_rates$p01
#> [1] 1
#> 
#> $weibull_rates$p02
#> [1] 1
#> 
#> $weibull_rates$p12
#> [1] 1
#> 
#> 
#> $family
#> [1] "exponential"
#> 
#> attr(,"class")
#> [1] "ExponentialTransition" "TransitionParameters" 
```
