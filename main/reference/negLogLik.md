# Compute the Negative Log-Likelihood for a Given Data Set and Transition Model

Compute the Negative Log-Likelihood for a Given Data Set and Transition
Model

## Usage

``` r
negLogLik(transition, data)
```

## Arguments

- transition:

  (`ExponentialTransition` or `WeibullTransition`)\
  see
  [`exponential_transition()`](https://insightsengineering.github.io/simIDM/reference/exponential_transition.md)
  or
  [`weibull_transition()`](https://insightsengineering.github.io/simIDM/reference/weibull_transition.md)
  for details.

- data:

  (`data.frame`)\
  in the format created by
  [`prepareData()`](https://insightsengineering.github.io/simIDM/reference/prepareData.md).

## Value

The value of the negative log-likelihood.

## Details

Calculates the negative log-likelihood for a given data set and
transition model. It uses the hazard and survival functions specific to
the transition model.

## Examples

``` r
transition <- exponential_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6)
simData <- getOneClinicalTrial(
  nPat = c(30), transitionByArm = list(transition),
  dropout = list(rate = 0.8, time = 12),
  accrual = list(param = "time", value = 1)
)
negLogLik(transition, prepareData(simData))
#> [1] 24.68872
```
