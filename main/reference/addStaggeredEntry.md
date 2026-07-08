# Staggered Study Entry

This function adds staggered study entry times to a simulated data set
with illness-death model structure.

## Usage

``` r
addStaggeredEntry(simData, N, accrualParam, accrualValue)
```

## Arguments

- simData:

  (`data.frame`)\
  simulated data frame containing entry and exit times at individual
  study time scale. See
  [`getSimulatedData()`](https://insightsengineering.github.io/simIDM/reference/getSimulatedData.md)
  for details.

- N:

  (`int`)\
  number of patients.

- accrualParam:

  (`string`)\
  possible values are 'time' or 'intensity'.

- accrualValue:

  (`number`)\
  specifies the accrual intensity. For `accrualParam` equal time, it
  describes the length of the accrual period. For `accrualParam` equal
  intensity, it describes the number of patients recruited per time
  unit. If `accrualValue` is equal to 0, all patients start at calendar
  time 0 in the initial state.

## Value

This returns a data set containing a single simulated study containing
accrual times, i.e. staggered study entry. This is a helper function of
[`getSimulatedData()`](https://insightsengineering.github.io/simIDM/reference/getSimulatedData.md).

## Examples

``` r
simData <- data.frame(
  id = c(1, 1, 2, 3), from = c(0, 1, 0, 0), to = c(1, 2, "cens", 2),
  entry = c(0, 3, 0, 0),
  exit = c(3, 5.3, 5.6, 7.2), censTime = c(6.8, 6.8, 5.6, 9.4)
)
addStaggeredEntry(simData, 3, accrualParam = "time", accrualValue = 5)
#>   id from   to entry exit  entryAct  exitAct   censAct
#> 1  1    0    1     0  3.0 4.1716652 7.171665 10.971665
#> 2  1    1    2     3  5.3 7.1716652 9.471665 10.971665
#> 3  2    0 cens     0  5.6 3.0038044 8.603804  8.603804
#> 4  3    0    2     0  7.2 0.7860422 7.986042 10.186042
```
