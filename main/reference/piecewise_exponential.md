# Transition Hazards for Piecewise Exponential Event Times

This creates a list with class `TransitionParameters` containing
hazards, time intervals and Weibull rates for piecewise exponential
event times in an illness-death model.

## Usage

``` r
piecewise_exponential(h01, h02, h12, pw01, pw02, pw12)
```

## Arguments

- h01:

  (`numeric vector`)\
  constant transition hazards for 0 to 1 transition

- h02:

  (`numeric vector`)\
  constant transition hazards for 0 to 2 transition

- h12:

  (`numeric vector`)\
  constant transition hazards for 1 to 2 transition

- pw01:

  (`numeric vector`)\
  time intervals for the piecewise constant hazards `h01`

- pw02:

  (`numeric vector`)\
  time intervals for the piecewise constant hazards `h02`

- pw12:

  (`numeric vector`)\
  time intervals for the piecewise constant hazards `h12`

## Value

List with elements `hazards`, `intervals`, `weibull_rates` and `family`
(piecewise exponential).

## Examples

``` r
piecewise_exponential(
  h01 = c(1, 1, 1), h02 = c(1.5, 0.5, 1), h12 = c(1, 1, 1),
  pw01 = c(0, 3, 8), pw02 = c(0, 6, 7), pw12 = c(0, 8, 9)
)
#> $hazards
#> $hazards$h01
#> [1] 1 1 1
#> 
#> $hazards$h02
#> [1] 1.5 0.5 1.0
#> 
#> $hazards$h12
#> [1] 1 1 1
#> 
#> 
#> $intervals
#> $intervals$pw01
#> [1] 0 3 8
#> 
#> $intervals$pw02
#> [1] 0 6 7
#> 
#> $intervals$pw12
#> [1] 0 8 9
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
#> [1] "piecewise exponential"
#> 
#> attr(,"class")
#> [1] "PWCTransition"        "TransitionParameters"
```
