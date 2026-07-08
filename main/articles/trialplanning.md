# Power and Type I Error Calculations

## Introduction

This vignette demonstrates how to use `simIDM` for trial design planning
with a simple example. We will show how to estimate type I errors and
statistical power from simulations to optimize study design, details can
be found in Erdmann et al. (2023). Jointly modeling the endpoints PFS
and OS with the illness-death model has two major advantages: - We
properly account for the correlation between PFS and OS, - The
assumption of proportional hazards is not required.

OS is defined as the time to reach the absorbing state death, and PFS is
defined as the time to reach the absorbing state `death` or the
intermediate state `progression`, whichever occurs first. Figure 1 shows
the multistate model with the corresponding transition hazards. In the
vignette, we show how to estimate type I errors and statistical power
from simulations and give an idea on how this can be used to plan
complex study trials.

![ Figure 1 - Multistate model with indermediate state progession and
absorbing state
death](trialplanning_files/figure-html/unnamed-chunk-1-1.png)

Figure 1 - Multistate model with indermediate state progession and
absorbing state death

## Scenario - PFS and OS as Co-primary Endpoints

We consider the following study design:

- PFS and OS as co-primary endpoints with one final analysis each,
  i.e. for a successful trial all endpoints need to be significant.
- Treatment vs. control group, 1:1 randomization ratio
- Global significance level of 5 %
- The standard log-rank test is used to assess the null hypothesis of
  equal survival functions in both groups for PFS and OS, respectively.
- Statistical power to detect a difference between the groups should be
  80 % for each endpoint
- 5 % drop-out within 12 time units
- Accrual of 100 patients per time unit

Using the multistate model approach implies that trial planning is based
on assumptions on the three transition hazards in each arm, i.e. six
hazards in total. In our example scenario, we assume that all six
transition hazards are constant and a small effect of the treatment on
hazards from the initial state to death. Transition hazards are tuned
such that median time until a PFS event is 0.87 time units in the
control arm and 1.2 time units in the treatment arm. Median time until
an OS event is 1.76 time units in the control group and 2.1 time units
in the treatment group.

Figure 2 shows the transition hazards, survival functions, hazard
functions and hazard ratios for both endpoints.

![Figure 2 - Transition hazards, survival functions, hazard functions
and hazard ratios for our scenario.](scenario.png)

Figure 2 - Transition hazards, survival functions, hazard functions and
hazard ratios for our scenario.

The transition hazards are specified as follows:

``` r

library(simIDM)
transitionTrt <- exponential_transition(h01 = 0.3, h02 = 0.28, h12 = 0.5)
transitionPl <- exponential_transition(h01 = 0.5, h02 = 0.3, h12 = 0.6)

transitionList <- list(transitionPl, transitionTrt)
```

The package provides functions that return the values of the PFS or OS
survival functions for given transition hazards (Constant, Weibull or
Piecewise Constant) and pre-specified time points.

``` r

timepoints <- c(0, 0.1, 0.3, 0.7, 1, 5)
# OS Survival function for Constant transition hazards:
ExpSurvOS(timepoints, h01 = 0.2, h02 = 0.4, h12 = 0.1)
#> [1] 1.0000000 0.9610787 0.8893403 0.7671856 0.6912219 0.2724845
# OS Survival function for Weibull transition hazards:
WeibSurvOS(timepoints, h01 = 0.2, h02 = 0.5, h12 = 2.1, p01 = 1.2, p02 = 0.9, p12 = 1)
#> [1] 1.00000000 0.93822237 0.83706585 0.66353708 0.55296799 0.03684786
# OS Survival function for Piecewise Constant transition hazards:
PWCsurvOS(timepoints,
  h01 = c(0.3, 0.5), h02 = c(0.5, 0.8), h12 = c(0.7, 1),
  pw01 = c(0, 4), pw02 = c(0, 8), pw12 = c(0, 3)
)
#> [1] 1.00000000 0.95094877 0.85849702 0.69546105 0.59109798 0.03945673
```

There are also functions for PFS survival functions available:

``` r

timepoints <- c(0, 0.1, 0.3, 0.7, 1, 5)
# PFS Survival function for Constant transition hazards:
ExpSurvPFS(timepoints, h01 = 0.2, h02 = 0.4)
#> [1] 1.00000000 0.94176453 0.83527021 0.65704682 0.54881164 0.04978707
# PFS Survival function for Weibull transition hazards:
WeibSurvPFS(timepoints, h01 = 0.2, h02 = 0.5, p01 = 1.2, p02 = 0.9)
#> [1] 1.00000000 0.92721907 0.80545180 0.61074857 0.49658530 0.02995439
# PFS Survival function for Piecewise Constant transition hazards:
PWCsurvPFS(timepoints,
  h01 = c(0.3, 0.5), h02 = c(0.5, 0.8),
  pw01 = c(0, 4), pw02 = c(0, 8)
)
#> [1] 1.00000000 0.92311635 0.78662786 0.57120906 0.44932896 0.01499558
```

For PFS, the hazard ratio under $`H_0`$ is known by specification:

``` r

hTrtPFS <- sum(unlist(transitionTrt$hazards[1:2]))
hPlPFS <- sum(unlist(transitionPl$hazards[1:2]))
hRatioPFS <- hTrtPFS / hPlPFS
hRatioPFS
#> [1] 0.725
```

For OS, the ratio of hazard functions is not necessarily constant. An
averaged HR can be calculated using `avgHRExpOS`:

``` r

hRatioOS <- avgHRExpOS(transitionByArm = transitionList, alpha = 0.5, upper = 1000)
hRatioOS
#> [1] 0.8072368
```

## Type I Error - Simulation Under $`H_0`$

The type I error can be estimated empirically by simulating clinical
trials under $`H_0`$. To achieve this, we set the transition hazards of
the treatment group to match those of the control group. Then, we use
[`getClinicalTrials()`](https://insightsengineering.github.io/simIDM/reference/getClinicalTrials.md)
to generate a large number of simulated trials. We will use 100
iterations here. For applications, to achieve satisfactory precision in
estimates of type I error, a higher number (e.g. 10,000) is recommended.

``` r

transitionListNull <- list(transitionPl, transitionPl)
nRep <- 100
simNull <- getClinicalTrials(
  nRep = nRep, nPat = c(800, 800), seed = 1238, datType = "1rowPatient",
  transitionByArm = transitionListNull,
  dropout = list(rate = 0.05, time = 12), accrual = list(param = "intensity", value = 100)
)
```

Using the simulation, we can now identify critical log-rank test values
for both PFS and OS to maintain a 5% global significance level.
Initially, we allocate critical values so that the two-sided log-rank
test has a 4% significance level for the OS endpoint and 1% for the PFS
endpoint, effectively splitting the global significance level. Such a
Bonferroni correction is widely used in trials with co-primary
endpoints.

``` r

alphaOS <- 0.04
alphaPFS <- 0.01
criticalOS <- abs(qnorm(alphaOS / 2))
criticalPFS <- abs(qnorm(alphaPFS / 2))
```

With the Schoenfeld approximation, preliminary sample sizes can be
computed to get an idea of how many events are needed to achieve 80 %
power:

``` r

library(rpact)
#> Installation qualification for rpact 4.4.0 has not yet been performed.
#> Please run testPackage() before using the package in GxP relevant environments.
# Number of PFS events required for 80 % power via Schoenfeld:
schoenfeldPFS <- getSampleSizeSurvival(
  lambda2 = hPlPFS, hazardRatio = hRatioPFS,
  accrualTime = 0, accrualIntensity = 500,
  maxNumberOfSubjects = 1000, sided = 1,
  alpha = alphaPFS / 2, beta = 0.2
)
eventNumPFS <- ceiling(schoenfeldPFS$eventsFixed)
eventNumPFS
#> [1] 452

# Number of OS events required for 80 % power via Schoenfeld:
schoenfeldOS <- getSampleSizeSurvival(
  lambda2 = hPlPFS, hazardRatio = hRatioOS,
  accrualTime = 0, accrualIntensity = 500,
  maxNumberOfSubjects = 1000, sided = 1,
  alpha = alphaOS / 2, beta = 0.2
)
eventNumOS <- ceiling(schoenfeldOS$eventsFixed)
eventNumOS
#> [1] 732
```

Using these critical values and required number of events for PFS and
OS, we can now determine the global type I error empirically by counting
the number of trials simulated under $`H_0`$ with significant log-rank
tests. The empirical type I error for each endpoint is calculated as the
proportion of trials with significant log-rank tests, while the global
Type I error is the proportion of trials where at least one of PFS or OS
with a significant log-rank test. This can be done using
`empSignificant`:

``` r

empSignificant(
  simTrials = simNull,
  criticalPFS = criticalPFS,
  criticalOS = criticalOS,
  eventNumPFS = eventNumPFS,
  eventNumOS = eventNumOS
)
#> $significantPFS
#> [1] 0
#> 
#> $significantOS
#> [1] 0.03
#> 
#> $significantAtLeastOne
#> [1] 0.03
#> 
#> $significantBoth
#> [1] 0
```

In this example, the global type I error rate is 3%.

## Sample size and power calculation - simulation under $`H_1`$

Next, we simulate a large number of trials under $`H_1`$ to compute the
empirical power:

``` r

simH1 <- getClinicalTrials(
  nRep = nRep, nPat = c(800, 800), seed = 1238, datType = "1rowPatient",
  transitionByArm = transitionList,
  dropout = list(rate = 0.05, time = 12), accrual = list(param = "intensity", value = 100)
)
```

The empirical power for each endpoint is the proportion of simulated
trials with significant log-rank tests under $`H_1`$. The multistate
model approach allows us to easily estimate further interesting metrics,
such as joint power, i.e. the probability that both endpoints in a trial
are significant, if each endpoint is analyzed at its planned time-point.

``` r

empSignificant(
  simTrials = simH1,
  criticalPFS = criticalPFS,
  criticalOS = criticalOS,
  eventNumPFS = eventNumPFS,
  eventNumOS = eventNumOS
)
#> $significantPFS
#> [1] 0.83
#> 
#> $significantOS
#> [1] 0.75
#> 
#> $significantAtLeastOne
#> [1] 0.9
#> 
#> $significantBoth
#> [1] 0.68
```

In this example, for the endpoint OS, the number of events has to be
increased to obtain a power of 80 %. Similarly, we can reduce the number
of events if the simulation shows that a trial design appears to be
overpowered.

It is also possible to derive the median time at which a certain number
of events are expected to occur and how many events of the second
endpoint have occurred at that time on average.

``` r

# median time point for 329 PFS events to have occurred:
timePointsPFS <- lapply(simH1, getTimePoint,
  eventNum = 329, typeEvent = "PFS",
  byArm = FALSE
)
median(unlist(timePointsPFS))
#> [1] 2.912601

# median time point for 684 OS events to have occurred:
timePointsOS <- lapply(simH1, getTimePoint,
  eventNum = 684, typeEvent = "OS",
  byArm = FALSE
)
median(unlist(timePointsOS))
#> [1] 5.804041

# mean number of PFS events at time of OS analysis
eventsPFS <- lapply(
  seq_along(timePointsPFS),
  function(t) {
    sum(simH1[[t]]$OSevent[(simH1[[t]]$OStime + simH1[[t]]$recruitTime) <= timePointsPFS[[t]]])
  }
)
mean(unlist(eventsPFS))
#> [1] 219.31

# mean number of OS events at time of PFS analysis
eventsOS <- lapply(
  seq_along(timePointsOS),
  function(t) {
    sum(simH1[[t]]$PFSevent[(simH1[[t]]$PFStime + simH1[[t]]$recruitTime) <= timePointsOS[[t]]])
  }
)
mean(unlist(eventsOS))
#> [1] 864.08
```

## References

Erdmann, Alexandra, Jan Beyersmann, and Kaspar Rufibach. 2023. “Oncology
Clinical Trial Design Planning Based on a Multistate Model That Jointly
Models Progression-Free and Overall Survival Endpoints.” *arXiv Preprint
arXiv:2301.10059*.
