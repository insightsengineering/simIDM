# Single Piecewise Exponentially Distributed Event Time

This returns an event time with a distribution resulting from piece-wise
constant hazards using the inversion method.

## Usage

``` r
PCWInversionMethod(haz, pw, LogU)
```

## Arguments

- haz:

  (`numeric`)\
  piecewise constant hazard.

- pw:

  (`numeric`)\
  time intervals for the piecewise constant hazard.

- LogU:

  (`numeric`)\
  transformed uniformly distributed random variables (log(1-U)).

## Value

This returns one single event time.

## Examples

``` r
PCWInversionMethod(haz = c(1.1, 0.5, 0.4), pw = c(0, 7, 10), LogU = log(1 - runif(1)))
#> [1] 0.07654301
```
