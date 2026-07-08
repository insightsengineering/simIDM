# Helper Function for `survPFSOS()`

Helper Function for
[`survPFSOS()`](https://insightsengineering.github.io/simIDM/reference/survPFSOS.md)

## Usage

``` r
PFSOSInteg(u, t, transition)
```

## Arguments

- u:

  (`numeric`)\
  variable of integration.

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

Numeric result of the integrand used to calculate the PFS\*OS survival
function.

## Note

Not all vectors `u` and `t` work here due to assertions in
[`log_p11()`](https://insightsengineering.github.io/simIDM/reference/log_p11.md).
