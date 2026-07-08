# Survival Function for Different Transition Models

Survival Function for Different Transition Models

## Usage

``` r
survTrans(transition, t, trans)

# S3 method for class 'ExponentialTransition'
survTrans(transition, t, trans)

# S3 method for class 'WeibullTransition'
survTrans(transition, t, trans)
```

## Arguments

- transition:

  (`ExponentialTransition` or `WeibullTransition`)\
  see
  [`exponential_transition()`](https://insightsengineering.github.io/simIDM/reference/exponential_transition.md)
  or
  [`weibull_transition()`](https://insightsengineering.github.io/simIDM/reference/weibull_transition.md)
  for details.

- t:

  (`numeric`)\
  time at which survival probability is to be computed.

- trans:

  (`integer`)\
  index specifying the transition type.

## Value

The survival probability for the specified transition and time.

## Methods (by class)

- `survTrans(ExponentialTransition)`: for the Exponential Transition
  Model

- `survTrans(WeibullTransition)`: for the Weibull Transition Model

## Examples

``` r
transition <- exponential_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6)
survTrans(transition, 0.4, 2)
#> [1] 0.5488116
transition <- exponential_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6)
survTrans(transition, 0.4, 2)
#> [1] 0.5488116
transition <- weibull_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6, p01 = 2, p02 = 2.5, p12 = 3)
survTrans(transition, 0.4, 2)
#> [1] 0.8591693
```
