# Simulate Data Set from an Illness-Death Model

This function creates a single simulated data set for a single treatment
arm. It simulates data from an illness-death model with one row per
transition and subject.

## Usage

``` r
getSimulatedData(
  N,
  transition = exponential_transition(h01 = 1, h02 = 1, h12 = 1),
  dropout = list(rate = 0, time = 12),
  accrual = list(param = "time", value = 0)
)
```

## Arguments

- N:

  (`int`)\
  number of patients.

- transition:

  (`TransitionParameters`)\
  transition parameters comprising `hazards`, corresponding `intervals`
  and `weibull_rates`, see
  [`exponential_transition()`](https://insightsengineering.github.io/simIDM/reference/exponential_transition.md),
  [`piecewise_exponential()`](https://insightsengineering.github.io/simIDM/reference/piecewise_exponential.md)
  and
  [`weibull_transition()`](https://insightsengineering.github.io/simIDM/reference/weibull_transition.md)
  for details.

- dropout:

  (`list`)\
  specifies drop-out probability. Random censoring times are generated
  using exponential distribution. `dropout$rate` specifies the drop-out
  probability per `dropout$time` time units. If `dropout$rate` is equal
  to 0, then no censoring is applied.

- accrual:

  (`list`)\
  specifies accrual intensity. See
  [`addStaggeredEntry()`](https://insightsengineering.github.io/simIDM/reference/addStaggeredEntry.md)
  for details.

## Value

This returns a data frame with one row per transition per individual.

## Details

The output data set contains the following columns:

- id (`integer`): patient id.

- from (`numeric`): starting state of the transition.

- to (`character`): final state of the transition.

- entry (`numeric`): entry time of the transition on the individual time
  scale.

- exit (`numeric`): exit time of the transition on the individual time
  scale.

- entryAct (`numeric`): entry time of the transition on study time
  scale.

- exitAct (`numeric`): exit time of the transition on study time scale.

- censAct (`numeric`): censoring time of the individual on study time
  scale.

## Examples

``` r
getSimulatedData(
  N = 10,
  transition = exponential_transition(h01 = 1, h02 = 1.5, h12 = 1),
  dropout = list(rate = 0.3, time = 1),
  accrual = list(param = "time", value = 5)
)
#>    id from to     entry        exit  entryAct   exitAct   censAct
#> 1   1    0  1 0.0000000 0.620109694 2.2084749 2.8285846  4.641584
#> 2   1    1  2 0.6201097 1.133229600 2.8285846 3.3417045  4.641584
#> 3   2    0  2 0.0000000 0.008524979 4.7083564 4.7168814  7.666176
#> 4   3    0  1 0.0000000 0.179999738 4.4569281 4.6369278 10.992647
#> 5   3    1  2 0.1799997 2.352246465 4.6369278 6.8091746 10.992647
#> 6   4    0  2 0.0000000 0.100762491 0.7384638 0.8392262  9.039251
#> 7   5    0  2 0.0000000 0.269376008 2.7693123 3.0386883  3.476619
#> 8   6    0  2 0.0000000 0.251356782 4.2223962 4.4737530  9.877679
#> 9   7    0  2 0.0000000 0.895647415 2.6340240 3.5296714  7.783116
#> 10  8    0  2 0.0000000 0.099131641 3.7638254 3.8629571  4.624288
#> 11  9    0  1 0.0000000 0.522817343 1.2393274 1.7621448  3.094538
#> 12  9    1  2 0.5228173 1.329566488 1.7621448 2.5688939  3.094538
#> 13 10    0  1 0.0000000 0.126395549 3.1352993 3.2616948 16.829769
#> 14 10    1  2 0.1263955 0.478119852 3.2616948 3.6134191 16.829769
```
