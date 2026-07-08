# PFS Survival Function for Different Transition Models

PFS Survival Function for Different Transition Models

## Usage

``` r
survPFS(transition, t)

# S3 method for class 'ExponentialTransition'
survPFS(transition, t)

# S3 method for class 'WeibullTransition'
survPFS(transition, t)

# S3 method for class 'PWCTransition'
survPFS(transition, t)
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

- t:

  (`numeric`)\
  time at which the value of the PFS survival function is to be
  computed.

## Value

The value of the survival function for the specified transition and
time.

## Methods (by class)

- `survPFS(ExponentialTransition)`: Survival Function for an exponential
  transition model.

- `survPFS(WeibullTransition)`: Survival Function for a Weibull
  transition model.

- `survPFS(PWCTransition)`: Survival Function for a piecewise constant
  transition model.

## Examples

``` r
transition <- exponential_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6)
survPFS(transition, 0.4)
#> [1] 0.3395955
transition <- exponential_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6)
survPFS(transition, 0.4)
#> [1] 0.3395955
transition <- weibull_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6, p01 = 2, p02 = 2.5, p12 = 3)
survPFS(transition, 0.4)
#> [1] 0.7090783
transition <- piecewise_exponential(
  h01 = c(1, 1, 1), h02 = c(1.5, 0.5, 1), h12 = c(1, 1, 1),
  pw01 = c(0, 3, 8), pw02 = c(0, 6, 7), pw12 = c(0, 8, 9)
)
survPFS(transition, 0.4)
#> [1] 0.3678794
```
