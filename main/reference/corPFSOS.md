# Correlation of PFS and OS event times for data from the IDM

Correlation of PFS and OS event times for data from the IDM

## Usage

``` r
corPFSOS(
  data,
  transition,
  bootstrap = TRUE,
  bootstrap_n = 100,
  conf_level = 0.95
)
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

- bootstrap:

  (`flag`)\
  if `TRUE` computes confidence interval via bootstrap.

- bootstrap_n:

  (`count`)\
  number of bootstrap samples.

- conf_level:

  (`proportion`)\
  confidence level for the confidence interval.

## Value

The correlation of PFS and OS.

## Examples

``` r
transition <- exponential_transition(h01 = 1.2, h02 = 1.5, h12 = 1.6)
data <- getClinicalTrials(
  nRep = 1, nPat = c(100), seed = 1234, datType = "1rowTransition",
  transitionByArm = list(transition), dropout = list(rate = 0.5, time = 12),
  accrual = list(param = "intensity", value = 7)
)[[1]]
#> Simulating 1 trials:
#>   |                                                                              |                                                                      |   0%  |                                                                              |======================================================================| 100%
corPFSOS(data, transition = exponential_transition(), bootstrap = FALSE)
#> $corPFSOS
#> [1] 0.6015421
#> 
if (FALSE) { # \dontrun{
corPFSOS(data, transition = exponential_transition(), bootstrap = TRUE)
} # }
```
