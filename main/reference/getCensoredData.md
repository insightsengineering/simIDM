# Helper function for `censoringByNumberEvents`

Helper function for `censoringByNumberEvents`

## Usage

``` r
getCensoredData(time, event, data)
```

## Arguments

- time:

  (`numeric`)\
  event times.

- event:

  (`numeric`)\
  event indicator.

- data:

  (`data.frame`)\
  data frame including patient id `id`, recruiting time `recruitTime`
  and individual censoring time `censTimeInd`.

## Value

This function returns a data frame with columns: event time, censoring
indicator, event indicator and event time in calendar time.

## Examples

``` r
transition1 <- weibull_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6, p01 = 0.8, p02 = 0.9, p12 = 1)
transition2 <- weibull_transition(h01 = 1, h02 = 1.3, h12 = 1.7, p01 = 1.1, p02 = 0.9, p12 = 1.1)

simStudy <- getOneClinicalTrial(
  nPat = c(20, 20), transitionByArm = list(transition1, transition2),
  dropout = list(rate = 0.3, time = 10),
  accrual = list(param = "time", value = 7)
)
simStudyWide <- getDatasetWideFormat(simStudy)
simStudyWide$censTimeInd <- 5 - simStudyWide$recruitTime
NotRecruited <- simStudyWide$id[simStudyWide$censTimeInd < 0]
censoredData <- simStudyWide[!(simStudyWide$id %in% NotRecruited), ]
getCensoredData(time = censoredData$OStime, event = censoredData$OSevent, data = censoredData)
#>           time Censored event   timeCal
#> 1  0.785039045        0     1 3.6470730
#> 2  0.264152292        0     1 1.1893622
#> 3  0.114419885        0     1 0.3878959
#> 4  1.000460067        0     1 3.8449832
#> 5  0.090251060        0     1 4.1901475
#> 6  0.057681178        0     1 0.5050084
#> 7  1.864437872        0     1 4.5753707
#> 8  0.505286431        0     1 4.0252529
#> 9  0.216705922        0     1 2.0222081
#> 10 0.021787680        0     1 3.0283585
#> 11 0.269505690        0     1 3.7295752
#> 12 1.979445101        0     1 4.4032408
#> 13 0.706481098        0     1 3.6007674
#> 14 0.154665230        0     1 4.3768401
#> 15 1.513935268        0     1 3.2915008
#> 16 0.653169302        0     1 1.6780308
#> 17 0.008262510        0     1 0.9143437
#> 18 0.810799371        0     1 3.2100054
#> 19 0.549038100        0     1 4.5680621
#> 20 0.407977908        0     1 3.6237299
#> 21 0.325864568        0     1 0.9782198
#> 22 0.036779895        0     1 4.2710633
#> 23 0.000675416        0     1 2.7206341
#> 24 1.276758059        0     1 3.1704302
#> 25 0.012196997        0     1 3.5368427
#> 26 0.929896017        0     1 1.0338115
#> 27 0.051693371        0     1 1.3494330
#> 28 0.012141099        0     1 4.7478093
#> 29 0.823576155        0     1 3.1440231
```
