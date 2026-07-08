# Quantile function for OS survival function induced by an illness-death model

Quantile function for OS survival function induced by an illness-death
model

## Usage

``` r
ExpQuantOS(q = 1/2, h01, h02, h12)
```

## Arguments

- q:

  (`numeric`)\
  quantile(s) at which to compute event time (q = 1 / 2 corresponds to
  median).

- h01:

  (`numeric vector`)\
  constant transition hazards for 0 to 1 transition.

- h02:

  (`numeric vector`)\
  constant transition hazards for 0 to 2 transition.

- h12:

  (`numeric vector`)\
  constant transition hazards for 1 to 2 transition.

## Value

This returns the time(s) t such that the OS survival function at t
equals q.

## Examples

``` r
ExpQuantOS(1 / 2, 0.2, 0.5, 2.1)
#> [1] 1.144539
```
