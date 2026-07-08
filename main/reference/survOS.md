# OS Survival Function for Different Transition Models

OS Survival Function for Different Transition Models

## Usage

``` r
survOS(transition, t)

# S3 method for class 'ExponentialTransition'
survOS(transition, t)

# S3 method for class 'WeibullTransition'
survOS(transition, t)

# S3 method for class 'PWCTransition'
survOS(transition, t)
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
  time at which the value of the OS survival function is to be computed.

## Value

The value of the survival function for the specified transition and
time.

## Methods (by class)

- `survOS(ExponentialTransition)`: Survival Function for an exponential
  transition model.

- `survOS(WeibullTransition)`: Survival Function for a Weibull
  transition model.

- `survOS(PWCTransition)`: Survival Function for a piecewise constant
  transition model.

## Examples

``` r
transition <- exponential_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6)
survOS(transition, 0.4)
#> [1] 0.5443558
transition <- exponential_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6)
survOS(transition, 0.4)
#> [1] 0.5443558
transition <- weibull_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6, p01 = 2, p02 = 2.5, p12 = 3)
survOS(transition, 0.4)
#> [1] 0.8627844
transition <- piecewise_exponential(
  h01 = c(1, 1, 1), h02 = c(1.5, 0.5, 1), h12 = c(1, 1, 1),
  pw01 = c(0, 3, 8), pw02 = c(0, 6, 7), pw12 = c(0, 8, 9)
)
survOS(transition, 0.4)
#> [1] 0.5695065
```
