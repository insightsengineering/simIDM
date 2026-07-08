# Number of recruited/censored/ongoing Patients.

Number of recruited/censored/ongoing Patients.

## Usage

``` r
getEventsAll(data, t)
```

## Arguments

- data:

  (`data.frame`)\
  illness-death data set in `1rowPatient` format.

- t:

  (`numeric`)\
  study time-point.

## Value

This function returns number of recruited patients, number of censored
and number of patients under observations.

## Examples

``` r
transition1 <- weibull_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6, p01 = 0.8, p02 = 0.9, p12 = 1)
transition2 <- weibull_transition(h01 = 1, h02 = 1.3, h12 = 1.7, p01 = 1.1, p02 = 0.9, p12 = 1.1)

simStudy <- getOneClinicalTrial(
  nPat = c(20, 20), transitionByArm = list(transition1, transition2),
  dropout = list(rate = 0.6, time = 10),
  accrual = list(param = "time", value = 0)
)
simStudyWide <- getDatasetWideFormat(simStudy)
getEventsAll(data = simStudyWide, t = 1.5)
#> Recruited  Censored  UnderObs 
#>        40         1         8 
```
