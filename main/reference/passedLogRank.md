# Helper function to conduct log-rank tests for either PFS or OS

This function evaluates the significance of either PFS or OS endpoints
for each trial in a list of trials, based on a pre-specified critical
value.

## Usage

``` r
passedLogRank(simTrials, typeEvent, eventNum, critical)
```

## Arguments

- simTrials:

  (`list`)\
  simulated trial data sets, see
  [`getClinicalTrials()`](https://insightsengineering.github.io/simIDM/reference/getClinicalTrials.md).

- typeEvent:

  (`string`)\
  endpoint to be evaluated, possible values are `PFS` and `OS`.

- eventNum:

  (`integer`)\
  number of events required to trigger analysis.

- critical:

  (positive `number`)\
  critical value of the log-rank test.

## Value

Logical vector indicating log-rank test significance for each trial.

## Examples

``` r
transition1 <- exponential_transition(h01 = 0.06, h02 = 0.3, h12 = 0.3)
transition2 <- exponential_transition(h01 = 0.1, h02 = 0.4, h12 = 0.3)
simTrials <- getClinicalTrials(
  nRep = 50, nPat = c(800, 800), seed = 1234, datType = "1rowPatient",
  transitionByArm = list(transition1, transition2), dropout = list(rate = 0.5, time = 12),
  accrual = list(param = "intensity", value = 7)
)
#> Simulating 50 trials:
#>   |                                                                              |                                                                      |   0%  |                                                                              |=                                                                     |   2%  |                                                                              |===                                                                   |   4%  |                                                                              |====                                                                  |   6%  |                                                                              |======                                                                |   8%  |                                                                              |=======                                                               |  10%  |                                                                              |========                                                              |  12%  |                                                                              |==========                                                            |  14%  |                                                                              |===========                                                           |  16%  |                                                                              |=============                                                         |  18%  |                                                                              |==============                                                        |  20%  |                                                                              |===============                                                       |  22%  |                                                                              |=================                                                     |  24%  |                                                                              |==================                                                    |  26%  |                                                                              |====================                                                  |  28%  |                                                                              |=====================                                                 |  30%  |                                                                              |======================                                                |  32%  |                                                                              |========================                                              |  34%  |                                                                              |=========================                                             |  36%  |                                                                              |===========================                                           |  38%  |                                                                              |============================                                          |  40%  |                                                                              |=============================                                         |  42%  |                                                                              |===============================                                       |  44%  |                                                                              |================================                                      |  46%  |                                                                              |==================================                                    |  48%  |                                                                              |===================================                                   |  50%  |                                                                              |====================================                                  |  52%  |                                                                              |======================================                                |  54%  |                                                                              |=======================================                               |  56%  |                                                                              |=========================================                             |  58%  |                                                                              |==========================================                            |  60%  |                                                                              |===========================================                           |  62%  |                                                                              |=============================================                         |  64%  |                                                                              |==============================================                        |  66%  |                                                                              |================================================                      |  68%  |                                                                              |=================================================                     |  70%  |                                                                              |==================================================                    |  72%  |                                                                              |====================================================                  |  74%  |                                                                              |=====================================================                 |  76%  |                                                                              |=======================================================               |  78%  |                                                                              |========================================================              |  80%  |                                                                              |=========================================================             |  82%  |                                                                              |===========================================================           |  84%  |                                                                              |============================================================          |  86%  |                                                                              |==============================================================        |  88%  |                                                                              |===============================================================       |  90%  |                                                                              |================================================================      |  92%  |                                                                              |==================================================================    |  94%  |                                                                              |===================================================================   |  96%  |                                                                              |===================================================================== |  98%  |                                                                              |======================================================================| 100%
passedLogRank(simTrials = simTrials, typeEvent = "PFS", eventNum = 300, critical = 2.4)
#>  [1]  TRUE  TRUE FALSE  TRUE  TRUE FALSE  TRUE FALSE  TRUE FALSE  TRUE  TRUE
#> [13] FALSE  TRUE  TRUE  TRUE  TRUE  TRUE FALSE  TRUE  TRUE FALSE  TRUE  TRUE
#> [25] FALSE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE FALSE FALSE  TRUE  TRUE
#> [37]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE FALSE FALSE FALSE  TRUE  TRUE
#> [49]  TRUE  TRUE
```
