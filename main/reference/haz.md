# Hazard Function for Different Transition Models

Hazard Function for Different Transition Models

## Usage

``` r
haz(transition, t, trans)

# S3 method for class 'ExponentialTransition'
haz(transition, t, trans)

# S3 method for class 'WeibullTransition'
haz(transition, t, trans)

# S3 method for class 'PWCTransition'
haz(transition, t, trans)
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
  time at which hazard is to be computed.

- trans:

  (`integer`)\
  index specifying the transition type.

## Value

The hazard rate for the specified transition and time.

## Details

The transition types are:

- `1`: Transition from state 0 (stable) to 1 (progression).

- `2`: Transition from state 0 (stable) to 2 (death).

- `3`: Transition from state 1 (progression) to 2 (death).

## Methods (by class)

- `haz(ExponentialTransition)`: for an exponential transition model.

- `haz(WeibullTransition)`: for the Weibull transition model.

- `haz(PWCTransition)`: for the piecewise constant transition model.

## Examples

``` r
transition <- exponential_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6)
haz(transition, 0.4, 2)
#> [1] 1.5
transition <- exponential_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6)
haz(transition, 0.4, 2)
#> [1] 1.5
transition <- weibull_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6, p01 = 2, p02 = 2.5, p12 = 3)
haz(transition, 0.4, 2)
#> [1] 0.9486833
transition <- piecewise_exponential(
  h01 = c(1, 1, 1), h02 = c(1.5, 0.5, 1), h12 = c(1, 1, 1),
  pw01 = c(0, 3, 8), pw02 = c(0, 6, 7), pw12 = c(0, 8, 9)
)
haz(transition, 6, 2)
#> [1] 0.5
```
