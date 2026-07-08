# Helper Function for `log_p11()`

Helper Function for
[`log_p11()`](https://insightsengineering.github.io/simIDM/reference/log_p11.md)

## Usage

``` r
p11Integ(x, transition)
```

## Arguments

- x:

  (`numeric`)\
  variable of integration.

- transition:

  (`TransitionParameters`)\
  see
  [`exponential_transition()`](https://insightsengineering.github.io/simIDM/reference/exponential_transition.md),
  [`weibull_transition()`](https://insightsengineering.github.io/simIDM/reference/weibull_transition.md)
  or
  [`piecewise_exponential()`](https://insightsengineering.github.io/simIDM/reference/piecewise_exponential.md)
  for details.

## Value

Hazard rate at the specified time for the transition from progression to
death.
