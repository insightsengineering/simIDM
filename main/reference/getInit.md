# Retrieve Initial Parameter Vectors for Likelihood Maximization

Retrieve Initial Parameter Vectors for Likelihood Maximization

## Usage

``` r
getInit(transition)

# S3 method for class 'ExponentialTransition'
getInit(transition)

# S3 method for class 'WeibullTransition'
getInit(transition)
```

## Arguments

- transition:

  (`ExponentialTransition` or `WeibullTransition`)\
  containing the initial parameters. See
  [`exponential_transition()`](https://insightsengineering.github.io/simIDM/reference/exponential_transition.md)
  or
  [`weibull_transition()`](https://insightsengineering.github.io/simIDM/reference/weibull_transition.md)
  for details.

## Value

The numeric vector of initial parameters for likelihood maximization.

## Methods (by class)

- `getInit(ExponentialTransition)`: for the Exponential Transition Model

- `getInit(WeibullTransition)`: for the Weibull Transition Model

## Examples

``` r
transition <- exponential_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6)
getInit(transition)
#> [1] 1.2 1.5 1.6
transition <- exponential_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6)
getInit(transition)
#> [1] 1.2 1.5 1.6
transition <- weibull_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6, p01 = 2, p02 = 2.5, p12 = 3)
getInit(transition)
#> [1] 1.2 1.5 1.6 2.0 2.5 3.0
```
