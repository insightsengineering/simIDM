# Probability of Remaining in Progression Between Two Time Points for Different Transition Models

Probability of Remaining in Progression Between Two Time Points for
Different Transition Models

## Usage

``` r
log_p11(transition, s, t)
```

## Arguments

- transition:

  (`TransitionParameters`)\
  see
  [`exponential_transition()`](https://insightsengineering.github.io/simIDM/reference/exponential_transition.md),
  [`weibull_transition()`](https://insightsengineering.github.io/simIDM/reference/weibull_transition.md)
  or
  [`piecewise_exponential()`](https://insightsengineering.github.io/simIDM/reference/piecewise_exponential.md)
  for details.

- s:

  (`numeric`)\
  lower time points.

- t:

  (`numeric`)\
  higher time points.

## Value

This returns the natural logarithm of the probability of remaining in
progression (state 1) between two time points, conditional on being in
state 1 at the lower time point.

## Examples

``` r
transition <- exponential_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6)
log_p11(transition, 1, 3)
#> [1] -3.2
```
