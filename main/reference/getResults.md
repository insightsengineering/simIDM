# Format Results of Parameter Estimation for Different Transition Models

Format Results of Parameter Estimation for Different Transition Models

## Usage

``` r
getResults(transition, res)

# S3 method for class 'ExponentialTransition'
getResults(transition, res)

# S3 method for class 'WeibullTransition'
getResults(transition, res)
```

## Arguments

- transition:

  (`TransitionParameters`)\
  see
  [`exponential_transition()`](https://insightsengineering.github.io/simIDM/reference/exponential_transition.md)
  or
  [`weibull_transition()`](https://insightsengineering.github.io/simIDM/reference/weibull_transition.md)
  for details.

- res:

  (`numeric` vector)\
  vector of parameter estimates from the likelihood maximization
  procedure.

## Value

Returns a `TransitionParameters` object with parameter estimates.

## Methods (by class)

- `getResults(ExponentialTransition)`: for the Exponential Transition
  Model

- `getResults(WeibullTransition)`: for the Weibull Transition Model

## Examples

``` r
results <- c(1.2, 1.5, 1.6)
getResults(exponential_transition(), results)
#> $hazards
#> $hazards$h01
#> [1] 1.2
#> 
#> $hazards$h02
#> [1] 1.5
#> 
#> $hazards$h12
#> [1] 1.6
#> 
#> 
#> $intervals
#> $intervals$pw01
#> [1] 0
#> 
#> $intervals$pw02
#> [1] 0
#> 
#> $intervals$pw12
#> [1] 0
#> 
#> 
#> $weibull_rates
#> $weibull_rates$p01
#> [1] 1
#> 
#> $weibull_rates$p02
#> [1] 1
#> 
#> $weibull_rates$p12
#> [1] 1
#> 
#> 
#> $family
#> [1] "exponential"
#> 
#> attr(,"class")
#> [1] "ExponentialTransition" "TransitionParameters" 
results <- c(1.2, 1.5, 1.6)
getResults(exponential_transition(), results)
#> $hazards
#> $hazards$h01
#> [1] 1.2
#> 
#> $hazards$h02
#> [1] 1.5
#> 
#> $hazards$h12
#> [1] 1.6
#> 
#> 
#> $intervals
#> $intervals$pw01
#> [1] 0
#> 
#> $intervals$pw02
#> [1] 0
#> 
#> $intervals$pw12
#> [1] 0
#> 
#> 
#> $weibull_rates
#> $weibull_rates$p01
#> [1] 1
#> 
#> $weibull_rates$p02
#> [1] 1
#> 
#> $weibull_rates$p12
#> [1] 1
#> 
#> 
#> $family
#> [1] "exponential"
#> 
#> attr(,"class")
#> [1] "ExponentialTransition" "TransitionParameters" 
results <- c(1.2, 1.5, 1.6, 2, 2.5, 1)
getResults(weibull_transition(), results)
#> $hazards
#> $hazards$h01
#> [1] 1.2
#> 
#> $hazards$h02
#> [1] 1.5
#> 
#> $hazards$h12
#> [1] 1.6
#> 
#> 
#> $weibull_rates
#> $weibull_rates$p01
#> [1] 2
#> 
#> $weibull_rates$p02
#> [1] 2.5
#> 
#> $weibull_rates$p12
#> [1] 1
#> 
#> 
#> $intervals
#> $intervals$pw01
#> [1] 0
#> 
#> $intervals$pw02
#> [1] 0
#> 
#> $intervals$pw12
#> [1] 0
#> 
#> 
#> $family
#> [1] "Weibull"
#> 
#> attr(,"class")
#> [1] "WeibullTransition"    "TransitionParameters"
```
