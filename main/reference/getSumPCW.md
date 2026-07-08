# Sum of Two Piecewise Constant Hazards

This returns the sum of two piecewise constant hazards per interval.

## Usage

``` r
getSumPCW(haz1, haz2, pw1, pw2)
```

## Arguments

- haz1:

  (`numeric`)\
  first summand (piecewise constant hazard).

- haz2:

  (`numeric`)\
  second summand (piecewise constant hazard).

- pw1:

  (`numeric`)\
  time intervals of first summand.

- pw2:

  (`numeric`)\
  time intervals of second summand.

## Value

List with elements `hazards` and `intervals` for the sum of two
piecewise constant hazards.

## Examples

``` r
getSumPCW(c(1.2, 0.3, 0.6), c(1.2, 0.7, 1), c(0, 8, 9), c(0, 1, 4))
#> $hazards
#> [1] 2.4 1.9 2.2 1.3 1.6
#> 
#> $intervals
#> [1] 0 1 4 8 9
#> 
```
