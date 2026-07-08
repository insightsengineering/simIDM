# Piecewise Exponentially Distributed Event Times

This returns event times with a distribution resulting from piece-wise
constant hazards using the inversion method.

## Usage

``` r
getPCWDistr(U, haz, pw, t_0)
```

## Arguments

- U:

  (`numeric`)\
  uniformly distributed random variables.

- haz:

  (`numeric`)\
  piecewise constant hazard.

- pw:

  (`numeric`)\
  time intervals for the piecewise constant hazard.

- t_0:

  (`numeric`)\
  the starting times.

## Value

This returns a vector with event times.

## Examples

``` r
getPCWDistr(U = runif(3), haz = c(1.1, 0.5, 0.4), pw = c(0, 7, 10), t_0 = c(0, 1, 4.2))
#> [1] 0.338923881 0.265666910 0.001900963
```
