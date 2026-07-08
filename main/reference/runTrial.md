# Helper Function for Adding Progress Bar to Trial Simulation

Helper Function for Adding Progress Bar to Trial Simulation

## Usage

``` r
runTrial(x, pb, ...)
```

## Arguments

- x:

  (`int`)\
  iteration index within lapply.

- ...:

  parameters transferred to
  [`getOneClinicalTrial()`](https://insightsengineering.github.io/simIDM/reference/getOneClinicalTrial.md),
  see
  [`getOneClinicalTrial()`](https://insightsengineering.github.io/simIDM/reference/getOneClinicalTrial.md)
  for details.

## Value

This returns the same as
[`getOneClinicalTrial()`](https://insightsengineering.github.io/simIDM/reference/getOneClinicalTrial.md),
but updates the progress bar.
