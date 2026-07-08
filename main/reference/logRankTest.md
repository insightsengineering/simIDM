# Log-Rank Test for a Single Trial

This function evaluates the significance of either PFS or OS endpoints
in a trial, based on a pre-specified critical value.

## Usage

``` r
logRankTest(data, typeEvent = c("PFS", "OS"), critical)
```

## Arguments

- data:

  (`data.frame`)\
  data frame containing entry and exit times of an illness-death model.
  See
  [`getSimulatedData()`](https://insightsengineering.github.io/simIDM/reference/getSimulatedData.md)
  for details.

- typeEvent:

  (`string`)\
  endpoint to be evaluated, possible values are `PFS` and `OS`.

- critical:

  (positive `number`)\
  critical value of the log-rank test.

## Value

Logical value indicating log-rank test significance.

## Examples

``` r
transition1 <- exponential_transition(h01 = 0.06, h02 = 0.3, h12 = 0.3)
transition2 <- exponential_transition(h01 = 0.1, h02 = 0.4, h12 = 0.3)
simTrial <- getClinicalTrials(
  nRep = 1, nPat = c(800, 800), seed = 1234, datType = "1rowPatient",
  transitionByArm = list(transition1, transition2), dropout = list(rate = 0.5, time = 12),
  accrual = list(param = "intensity", value = 7)
)[[1]]
#> Simulating 1 trials:
#>   |                                                                              |                                                                      |   0%  |                                                                              |======================================================================| 100%
logRankTest(data = simTrial, typeEvent = "OS", critical = 3.4)
#> [1] TRUE
```
