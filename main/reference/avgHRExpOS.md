# Average OS Hazard Ratio from Constant Transition Hazards

Average OS Hazard Ratio from Constant Transition Hazards

## Usage

``` r
avgHRExpOS(transitionByArm, alpha = 0.5, upper = Inf)
```

## Arguments

- transitionByArm:

  (`list`)\
  transition parameters for each treatment group. See
  [`exponential_transition()`](https://insightsengineering.github.io/simIDM/reference/exponential_transition.md),
  [`piecewise_exponential()`](https://insightsengineering.github.io/simIDM/reference/piecewise_exponential.md)
  and
  [`weibull_transition()`](https://insightsengineering.github.io/simIDM/reference/weibull_transition.md)
  for details.

- alpha:

  (`number`)\
  assigns weights to time points, where values higher than 0.5 assign
  greater weight to earlier times and values lower than 0.5 assign
  greater weight to later times.

- upper:

  (`number`)\
  upper (time) limit of integration.

## Value

This returns the value of the average hazard ratio.

## Examples

``` r
transitionTrt <- exponential_transition(h01 = 0.18, h02 = 0.06, h12 = 0.17)
transitionCtl <- exponential_transition(h01 = 0.23, h02 = 0.07, h12 = 0.19)
transitionList <- list(transitionCtl, transitionTrt)
avgHRExpOS(transitionByArm = transitionList, alpha = 0.5, upper = 100)
#> [1] 0.8319853
```
