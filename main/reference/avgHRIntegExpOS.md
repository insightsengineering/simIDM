# Helper Function for `avgHRExpOS()`

It is an integrand of the form OS hazard function with intensities h01,
h02, h12 at time point t multiplied with a weighted product of the two
OS Survival functions at t (one for intensities h0 and one for h1).

## Usage

``` r
avgHRIntegExpOS(x, h01, h02, h12, h0, h1, alpha)
```

## Arguments

- x:

  (`numeric`)\
  variable of integration.

- h01:

  (positive `number`)\
  transition hazard for 0 to 1 transition.

- h02:

  (positive `number`)\
  transition hazard for 0 to 2 transition.

- h12:

  (positive `number`)\
  transition hazard for 1 to 2 transition.

- h0:

  (`list`)\
  transition parameters for the first treatment group.

- h1:

  (`list`)\
  transition parameters for the second treatment group.

- alpha:

  (`number`)\
  weight parameter, see
  [`avgHRExpOS()`](https://insightsengineering.github.io/simIDM/reference/avgHRExpOS.md).

## Value

This returns the value of the integrand used to calculate the average
hazard ratio for constant transition hazards, see
[`avgHRExpOS()`](https://insightsengineering.github.io/simIDM/reference/avgHRExpOS.md).

## Examples

``` r
h0 <- list(h01 = 0.18, h02 = 0.06, h12 = 0.17)
h1 <- list(h01 = 0.23, h02 = 0.07, h12 = 0.19)
avgHRIntegExpOS(x = 5, h01 = 0.2, h02 = 0.5, h12 = 0.7, h0 = h0, h1 = h1, alpha = 0.5)
#> [1] 0.297362
```
