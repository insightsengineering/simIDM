# Helper Function for Computing E(OS^2)

Helper Function for Computing E(OS^2)

## Usage

``` r
expvalOSInteg(x, transition)
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

Numeric results of the integrand used to calculate E(OS^2).

## Examples

``` r
transition <- exponential_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6)
expvalOSInteg(0.4, transition)
#> [1] 0.2177423
```
