# Event-driven censoring.

This function censors a study after a pre-specified number of events
occurred.

## Usage

``` r
censoringByNumberEvents(data, eventNum, typeEvent)
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

## Value

This function returns a data set that is censored after `eventNum` of
`typeEvent`-events occurred.

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
censoringByNumberEvents(data = simStudyWide, eventNum = 20, typeEvent = "PFS")
#>    id trt    PFStime PFSevent     OStime CensoredOS OSevent recruitTime
#> 1   1   1 0.20251761        1 0.20251761          0       1   2.0303511
#> 2   2   1 0.32225015        1 0.33107435          1       0   3.3605262
#> 3   4   1 0.02100692        1 0.02100692          0       1   2.8050413
#> 4   5   1 0.48344054        1 0.48344054          0       1   1.4922090
#> 5   7   1 1.90043827        1 1.90043827          0       1   0.4102988
#> 6   9   1 0.13622285        1 0.13622285          0       1   1.0432483
#> 7  10   1 0.06170412        0 0.06170412          1       0   3.6298964
#> 8  13   1 0.22595459        1 0.22595459          0       1   1.6891981
#> 9  16   1 0.05506586        1 0.05506586          0       1   0.1956922
#> 10 17   1 0.40591047        1 0.40591047          0       1   3.2856901
#> 11 20   1 0.59979587        1 0.59979587          0       1   2.8273770
#> 12 22   2 0.11168171        1 0.11168171          0       1   0.7045952
#> 13 25   2 0.51764752        1 0.67757502          0       1   2.6835947
#> 14 27   2 0.22308854        0 0.22308854          1       0   2.4450933
#> 15 29   2 0.25793778        1 0.25793778          0       1   1.5126998
#> 16 30   2 0.16932117        1 0.28035167          0       1   0.2246490
#> 17 31   2 0.35858824        1 1.98005050          0       1   1.0172109
#> 18 33   2 0.26420629        1 0.26420629          0       1   1.4920452
#> 19 34   2 1.08675230        1 1.29091031          0       1   1.4721752
#> 20 35   2 1.57561614        1 1.57561614          0       1   0.2766448
#> 21 37   2 1.01421133        1 1.25121844          0       1   1.7144959
#> 22 39   2 0.64035375        1 1.67394039          1       0   2.0176602
#>    OStimeCal PFStimeCal
#> 1  2.2328687  2.2328687
#> 2  3.6916006  3.6827764
#> 3  2.8260482  2.8260482
#> 4  1.9756495  1.9756495
#> 5  2.3107371  2.3107371
#> 6  1.1794711  1.1794711
#> 7  3.6916006  3.6916006
#> 8  1.9151527  1.9151527
#> 9  0.2507580  0.2507580
#> 10 3.6916006  3.6916006
#> 11 3.4271729  3.4271729
#> 12 0.8162769  0.8162769
#> 13 3.3611697  3.2012422
#> 14 2.6681819  2.6681819
#> 15 1.7706376  1.7706376
#> 16 0.5050006  0.3939701
#> 17 2.9972614  1.3757991
#> 18 1.7562515  1.7562515
#> 19 2.7630855  2.5589275
#> 20 1.8522610  1.8522610
#> 21 2.9657144  2.7287072
#> 22 3.6916006  2.6580139
```
