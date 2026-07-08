# Time-point by which a specified number of events occurred.

This returns the study time-point by which a specified number of events
(PFS or OS) occurred.

## Usage

``` r
getTimePoint(data, eventNum, typeEvent, byArm = FALSE)
```

## Arguments

- data:

  (`data.frame`)\
  illness-death data set in `1rowPatient` format.

- eventNum:

  (`int`)\
  number of events.

- typeEvent:

  (`string`)\
  type of event. Possible values are `PFS` and `OS`.

- byArm:

  (`logical`)\
  if `TRUE` time-point per treatment arm, else joint evaluation of
  treatment arms.

## Value

This returns the time-point by which `eventNum` of `typeEvent`-events
occurred.

## Examples

``` r
transition1 <- weibull_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6, p01 = 0.8, p02 = 0.9, p12 = 1)
transition2 <- weibull_transition(h01 = 1, h02 = 1.3, h12 = 1.7, p01 = 1.1, p02 = 0.9, p12 = 1.1)

simStudy <- getOneClinicalTrial(
  nPat = c(20, 20), transitionByArm = list(transition1, transition2),
  dropout = list(rate = 0.3, time = 10),
  accrual = list(param = "time", value = 0)
)
simStudyWide <- getDatasetWideFormat(simStudy)
getTimePoint(simStudyWide, eventNum = 10, typeEvent = "OS", byArm = FALSE)
#>      all 
#> 0.124277 
```
