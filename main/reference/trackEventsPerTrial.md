# Event tracking in an oncology trial.

Event tracking in an oncology trial.

## Usage

``` r
trackEventsPerTrial(data, timeP, byArm = FALSE)
```

## Arguments

- data:

  (`data.frame`)\
  illness-death data set in `1rowPatient` format.

- timeP:

  (`numeric`)\
  vector of study time-points.

- byArm:

  (`logical`)\
  if `TRUE` time-point per treatment arm, else joint evaluation of
  treatment arms.

## Value

This function returns a data frame including number of PFS events,
number of OS events, number of recruited patients, number of censored
patients and number of ongoing patients at `timeP`.

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
trackEventsPerTrial(data = simStudyWide, timeP = 1.5, byArm = FALSE)
#> $all
#>           Timepoint: 1.5
#> PFS                   38
#> OS                    36
#> Recruited             40
#> Censored               1
#> Ongoing                3
#> 
```
