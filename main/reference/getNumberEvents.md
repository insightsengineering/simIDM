# Helper Function for `trackEventsPerTrial`

Helper Function for `trackEventsPerTrial`

## Usage

``` r
getNumberEvents(event, time, t)
```

## Arguments

- event:

  (`numeric`)\
  event indicator.

- time:

  (`numeric`)\
  event times.

- t:

  (`numeric`)\
  study time-point.

## Value

This function returns the number of events occurred until time t.

## Examples

``` r
event <- c(0, 1, 1, 1, 0)
time <- c(3, 3.4, 5, 6, 5.5)
getNumberEvents(event = event, time = time, t = 5)
#> [1] 2
```
