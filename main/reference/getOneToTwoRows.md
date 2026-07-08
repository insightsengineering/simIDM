# Transitions from the Intermediate State to the Absorbing State

This function creates transition entry and exit times from the
intermediate state to the absorbing state for an existing data frame
containing the exit times out of the initial state.

## Usage

``` r
getOneToTwoRows(simDataOne, transition)
```

## Arguments

- simDataOne:

  (`data.frame`)\
  a data frame containing all patients with transitions into the
  intermediate state. See
  [`getSimulatedData()`](https://insightsengineering.github.io/simIDM/reference/getSimulatedData.md)
  for details.

- transition:

  (`TransitionParameters`)\
  transition parameters comprising `hazards`, corresponding `intervals`
  and `weibull_rates`, see
  [`exponential_transition()`](https://insightsengineering.github.io/simIDM/reference/exponential_transition.md),
  [`piecewise_exponential()`](https://insightsengineering.github.io/simIDM/reference/piecewise_exponential.md)
  and
  [`weibull_transition()`](https://insightsengineering.github.io/simIDM/reference/weibull_transition.md)
  for details.

## Value

This returns a data frame with one row per patient for the second
transition, i.e. the transition out of the intermediate state. This is a
helper function of
[`getSimulatedData()`](https://insightsengineering.github.io/simIDM/reference/getSimulatedData.md).

## Examples

``` r
simDataOne <- data.frame(
  id = c(1:3), to = c(1, 1, 1), from = c(0, 0, 0), entry = c(0, 0, 0),
  exit = c(3, 5.6, 7.2), censTime = c(6.8, 5.9, 9.4)
)
transition <- exponential_transition(1, 1.6, 0.3)
getOneToTwoRows(simDataOne, transition)
#>   id from   to entry     exit censTime
#> 1  1    1    2   3.0 4.001717      6.8
#> 2  2    1 cens   5.6 5.900000      5.9
#> 3  3    1    2   7.2 8.366805      9.4
```
