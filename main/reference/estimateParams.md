# Estimate Parameters of the Multistate Model Using Clinical Trial Data

Estimate Parameters of the Multistate Model Using Clinical Trial Data

## Usage

``` r
estimateParams(data, transition)
```

## Arguments

- data:

  (`data.frame`)\
  in the format produced by
  [`getOneClinicalTrial()`](https://insightsengineering.github.io/simIDM/reference/getOneClinicalTrial.md).

- transition:

  (`TransitionParameters` object)\
  specifying the assumed distribution of transition hazards. Initial
  parameters for optimization can be specified here. See
  [`exponential_transition()`](https://insightsengineering.github.io/simIDM/reference/exponential_transition.md)
  or
  [`weibull_transition()`](https://insightsengineering.github.io/simIDM/reference/weibull_transition.md)
  for details.

## Value

Returns a `TransitionParameters` object with the estimated parameters.

## Details

This function estimates parameters for transition models using clinical
trial data. The `transition` object can be initialized with starting
values for parameter estimation. It uses
[`stats::optim()`](https://rdrr.io/r/stats/optim.html) to optimize the
parameters.

## Examples

``` r
transition <- exponential_transition(h01 = 2, h02 = 1.4, h12 = 1.6)
simData <- getOneClinicalTrial(
  nPat = c(30), transitionByArm = list(transition),
  dropout = list(rate = 0.3, time = 12),
  accrual = list(param = "time", value = 1)
)
# Initialize transition with desired starting values for optimization:
transition <- exponential_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6)
estimate <- estimateParams(simData, transition)
```
