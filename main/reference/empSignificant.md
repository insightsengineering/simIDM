# Empirical Significance for a List of Simulated Trials

This function computes four types of empirical significance — PFS, OS,
at-least (significant in at least one of PFS/OS), and joint (significant
in both PFS and OS) — using the log-rank test. Empirical significance is
calculated as the proportion of significant results in simulated trials,
each ending when a set number of PFS/OS events occur. Critical values
for PFS and OS test significance must be specified. If trials simulate
equal transition hazards across groups (H0), empirical significance
estimates type I error; if they simulate differing transition hazards
(H1), it estimates power.

## Usage

``` r
empSignificant(simTrials, criticalPFS, criticalOS, eventNumPFS, eventNumOS)
```

## Arguments

- simTrials:

  (`list`)\
  simulated trial data sets, see
  [`getClinicalTrials()`](https://insightsengineering.github.io/simIDM/reference/getClinicalTrials.md).

- criticalPFS:

  (positive `number`)\
  critical value of the log-rank test for PFS.

- criticalOS:

  (positive `number`)\
  critical value of the log-rank test for OS.

- eventNumPFS:

  (`integer`)\
  number of PFS events required to trigger PFS analysis.

- eventNumOS:

  (`integer`)\
  number of OS events required to trigger OS analysis.

## Value

This returns values of four measures of empirical significance.

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
empSignificant(
  simTrials = simTrials, criticalPFS = 2.4, criticalOS = 2.2,
  eventNumPFS = 300, eventNumOS = 500
)
#> $significantPFS
#> [1] 0.74
#> 
#> $significantOS
#> [1] 0.52
#> 
#> $significantAtLeastOne
#> [1] 0.78
#> 
#> $significantBoth
#> [1] 0.48
#> 
```
