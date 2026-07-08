# Generate the Target Function for Optimization

Generate the Target Function for Optimization

## Usage

``` r
getTarget(transition)

# S3 method for class 'ExponentialTransition'
getTarget(transition)

# S3 method for class 'WeibullTransition'
getTarget(transition)
```

## Arguments

- transition:

  (`TransitionParameters`)\
  specifying the distribution family. See
  [`exponential_transition()`](https://insightsengineering.github.io/simIDM/reference/exponential_transition.md)
  or
  [`weibull_transition()`](https://insightsengineering.github.io/simIDM/reference/weibull_transition.md)
  for details.

## Value

Function that calculates the negative log-likelihood for the given
parameters.

## Details

This function creates a target function for optimization, computing the
negative log-likelihood for given parameters, data, and transition model
type.

## Methods (by class)

- `getTarget(ExponentialTransition)`: for the Exponential Transition
  Model

- `getTarget(WeibullTransition)`: for the Weibull Transition Model

## Examples

``` r
transition <- exponential_transition(2, 1.3, 0.8)
simData <- getOneClinicalTrial(
  nPat = c(30), transitionByArm = list(transition),
  dropout = list(rate = 0.8, time = 12),
  accrual = list(param = "time", value = 1)
)
params <- c(1.2, 1.5, 1.6) # For ExponentialTransition
data <- prepareData(simData)
transition <- exponential_transition()
fun <- getTarget(transition)
fun(params, data)
#> [1] 45.08181
transition <- exponential_transition(2, 1.3, 0.8)
simData <- getOneClinicalTrial(
  nPat = c(30), transitionByArm = list(transition),
  dropout = list(rate = 0.8, time = 12),
  accrual = list(param = "time", value = 1)
)
params <- c(1.2, 1.5, 1.6)
data <- prepareData(simData)
transition <- exponential_transition()
target <- getTarget(transition)
target(params, data)
#> [1] 48.55708
transition <- weibull_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6, p01 = 2, p02 = 2.5, p12 = 3)
simData <- getOneClinicalTrial(
  nPat = c(30), transitionByArm = list(transition),
  dropout = list(rate = 0.8, time = 12),
  accrual = list(param = "time", value = 1)
)
params <- c(1.2, 1.5, 1.6, 0.8, 1.3, 1.1)
data <- prepareData(simData)
transition <- weibull_transition()
target <- getTarget(transition)
target(params, data)
#> [1] 34.37237
```
