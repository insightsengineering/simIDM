# Helper Function for Single Quantile for OS Survival Function

Helper Function for Single Quantile for OS Survival Function

## Usage

``` r
singleExpQuantOS(q, h01, h02, h12)
```

## Arguments

- q:

  (`number`)\
  single quantile at which to compute event time.

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

Single time t such that the OS survival function at t equals q.
