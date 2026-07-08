# Correlation of PFS and OS event times for Different Transition Models

Correlation of PFS and OS event times for Different Transition Models

## Usage

``` r
corTrans(transition)
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

## Value

The correlation of PFS and OS.

## Examples

``` r
transition <- exponential_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6)
corTrans(transition)
#> [1] 0.580381
```
