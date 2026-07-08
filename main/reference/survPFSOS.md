# Survival Function of the Product PFS\*OS for Different Transition Models

Survival Function of the Product PFS\*OS for Different Transition Models

## Usage

``` r
survPFSOS(t, transition)
```

## Arguments

- t:

  (`numeric`)\
  time at which the value of the PFS\*OS survival function is to be
  computed.

- transition:

  (`TransitionParameters`)\
  see
  [`exponential_transition()`](https://insightsengineering.github.io/simIDM/reference/exponential_transition.md),
  [`weibull_transition()`](https://insightsengineering.github.io/simIDM/reference/weibull_transition.md)
  or
  [`piecewise_exponential()`](https://insightsengineering.github.io/simIDM/reference/piecewise_exponential.md)
  for details.

## Value

This returns the value of PFS\*OS survival function at time t.

## Examples

``` r
transition <- exponential_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6)
survPFSOS(0.4, transition)
#> [1] 0.2480364
```
