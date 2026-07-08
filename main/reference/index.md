# Package index

## Package

- [`simIDM`](https://insightsengineering.github.io/simIDM/reference/simIDM-package.md)
  [`simIDM-package`](https://insightsengineering.github.io/simIDM/reference/simIDM-package.md)
  :

  `simIDM` Package

## Assertions

- [`assert_positive_number()`](https://insightsengineering.github.io/simIDM/reference/assert_positive_number.md)
  : Assertion for Positive Number
- [`assert_intervals()`](https://insightsengineering.github.io/simIDM/reference/assert_intervals.md)
  : Assertion for vector describing intervals

## Transition parameters

- [`exponential_transition()`](https://insightsengineering.github.io/simIDM/reference/exponential_transition.md)
  : Transition Hazards for Exponential Event Times
- [`piecewise_exponential()`](https://insightsengineering.github.io/simIDM/reference/piecewise_exponential.md)
  : Transition Hazards for Piecewise Exponential Event Times
- [`weibull_transition()`](https://insightsengineering.github.io/simIDM/reference/weibull_transition.md)
  : Transition Hazards for Weibull Distributed Event Times

## Helper functions for piecewise or Weibull Hazards

- [`getPWCHazard()`](https://insightsengineering.github.io/simIDM/reference/getPWCHazard.md)
  : Piecewise Constant Hazard Values
- [`getSumPCW()`](https://insightsengineering.github.io/simIDM/reference/getSumPCW.md)
  : Sum of Two Piecewise Constant Hazards
- [`getPCWDistr()`](https://insightsengineering.github.io/simIDM/reference/getPCWDistr.md)
  : Piecewise Exponentially Distributed Event Times
- [`PCWInversionMethod()`](https://insightsengineering.github.io/simIDM/reference/PCWInversionMethod.md)
  : Single Piecewise Exponentially Distributed Event Time
- [`getWaitTimeSum()`](https://insightsengineering.github.io/simIDM/reference/getWaitTimeSum.md)
  : Event Times Distributed as Sum of Weibull

## Simulation of a single data set for a single treatment arm

- [`getSimulatedData()`](https://insightsengineering.github.io/simIDM/reference/getSimulatedData.md)
  : Simulate Data Set from an Illness-Death Model
- [`getOneToTwoRows()`](https://insightsengineering.github.io/simIDM/reference/getOneToTwoRows.md)
  : Transitions from the Intermediate State to the Absorbing State
- [`addStaggeredEntry()`](https://insightsengineering.github.io/simIDM/reference/addStaggeredEntry.md)
  : Staggered Study Entry

## Simulation of a large number of clinical trials

- [`getOneClinicalTrial()`](https://insightsengineering.github.io/simIDM/reference/getOneClinicalTrial.md)
  : Simulation of a Single Oncology Clinical Trial
- [`getDatasetWideFormat()`](https://insightsengineering.github.io/simIDM/reference/getDatasetWideFormat.md)
  : Conversion of a Data Set from One Row per Transition to One Row per
  Patient
- [`getClinicalTrials()`](https://insightsengineering.github.io/simIDM/reference/getClinicalTrials.md)
  : Simulation of a Large Number of Oncology Clinical Trials

## Event tracking

- [`getTimePoint()`](https://insightsengineering.github.io/simIDM/reference/getTimePoint.md)
  : Time-point by which a specified number of events occurred.

- [`getCensoredData()`](https://insightsengineering.github.io/simIDM/reference/getCensoredData.md)
  :

  Helper function for `censoringByNumberEvents`

- [`censoringByNumberEvents()`](https://insightsengineering.github.io/simIDM/reference/censoringByNumberEvents.md)
  : Event-driven censoring.

- [`getEventsAll()`](https://insightsengineering.github.io/simIDM/reference/getEventsAll.md)
  : Number of recruited/censored/ongoing Patients.

- [`getNumberEvents()`](https://insightsengineering.github.io/simIDM/reference/getNumberEvents.md)
  :

  Helper Function for `trackEventsPerTrial`

- [`trackEventsPerTrial()`](https://insightsengineering.github.io/simIDM/reference/trackEventsPerTrial.md)
  : Event tracking in an oncology trial.

## Survival functions

- [`ExpSurvPFS()`](https://insightsengineering.github.io/simIDM/reference/ExpSurvPFS.md)
  : PFS Survival Function from Constant Transition Hazards
- [`ExpSurvOS()`](https://insightsengineering.github.io/simIDM/reference/ExpSurvOS.md)
  : OS Survival Function from Constant Transition Hazards
- [`ExpQuantOS()`](https://insightsengineering.github.io/simIDM/reference/ExpQuantOS.md)
  : Quantile function for OS survival function induced by an
  illness-death model
- [`WeibSurvPFS()`](https://insightsengineering.github.io/simIDM/reference/WeibSurvPFS.md)
  : PFS Survival Function from Weibull Transition Hazards
- [`WeibSurvOS()`](https://insightsengineering.github.io/simIDM/reference/WeibSurvOS.md)
  : OS Survival Function from Weibull Transition Hazards
- [`pwA()`](https://insightsengineering.github.io/simIDM/reference/pwA.md)
  : Cumulative Hazard for Piecewise Constant Hazards
- [`PWCsurvPFS()`](https://insightsengineering.github.io/simIDM/reference/PWCsurvPFS.md)
  : PFS Survival Function from Piecewise Constant Hazards
- [`PWCsurvOS()`](https://insightsengineering.github.io/simIDM/reference/PWCsurvOS.md)
  : OS Survival Function from Piecewise Constant Hazards

## Hazard functions

- [`ExpHazOS()`](https://insightsengineering.github.io/simIDM/reference/ExpHazOS.md)
  : OS Hazard Function from Constant Transition Hazards

- [`avgHRIntegExpOS()`](https://insightsengineering.github.io/simIDM/reference/avgHRIntegExpOS.md)
  :

  Helper Function for
  [`avgHRExpOS()`](https://insightsengineering.github.io/simIDM/reference/avgHRExpOS.md)

- [`avgHRExpOS()`](https://insightsengineering.github.io/simIDM/reference/avgHRExpOS.md)
  : Average OS Hazard Ratio from Constant Transition Hazards

## Empirical significance

- [`logRankTest()`](https://insightsengineering.github.io/simIDM/reference/logRankTest.md)
  : Log-Rank Test for a Single Trial
- [`empSignificant()`](https://insightsengineering.github.io/simIDM/reference/empSignificant.md)
  : Empirical Significance for a List of Simulated Trials

## Parameter Estimation

- [`estimateParams()`](https://insightsengineering.github.io/simIDM/reference/estimateParams.md)
  : Estimate Parameters of the Multistate Model Using Clinical Trial
  Data
- [`getInit()`](https://insightsengineering.github.io/simIDM/reference/getInit.md)
  : Retrieve Initial Parameter Vectors for Likelihood Maximization
- [`getResults()`](https://insightsengineering.github.io/simIDM/reference/getResults.md)
  : Format Results of Parameter Estimation for Different Transition
  Models
- [`getTarget()`](https://insightsengineering.github.io/simIDM/reference/getTarget.md)
  : Generate the Target Function for Optimization
- [`prepareData()`](https://insightsengineering.github.io/simIDM/reference/prepareData.md)
  : Preparation of a Data Set to Compute Log-likelihood

## Likelihood

- [`haz()`](https://insightsengineering.github.io/simIDM/reference/haz.md)
  : Hazard Function for Different Transition Models
- [`survTrans()`](https://insightsengineering.github.io/simIDM/reference/survTrans.md)
  : Survival Function for Different Transition Models
- [`negLogLik()`](https://insightsengineering.github.io/simIDM/reference/negLogLik.md)
  : Compute the Negative Log-Likelihood for a Given Data Set and
  Transition Model

## PFS-OS Correlation

- [`corPFSOS()`](https://insightsengineering.github.io/simIDM/reference/corPFSOS.md)
  : Correlation of PFS and OS event times for data from the IDM
- [`corTrans()`](https://insightsengineering.github.io/simIDM/reference/corTrans.md)
  : Correlation of PFS and OS event times for Different Transition
  Models
- [`expvalPFSInteg()`](https://insightsengineering.github.io/simIDM/reference/expvalPFSInteg.md)
  : Helper Function for Computing E(PFS^2)
- [`expvalOSInteg()`](https://insightsengineering.github.io/simIDM/reference/expvalOSInteg.md)
  : Helper Function for Computing E(OS^2)
- [`log_p11()`](https://insightsengineering.github.io/simIDM/reference/log_p11.md)
  : Probability of Remaining in Progression Between Two Time Points for
  Different Transition Models
- [`survPFS()`](https://insightsengineering.github.io/simIDM/reference/survPFS.md)
  : PFS Survival Function for Different Transition Models
- [`survPFSOS()`](https://insightsengineering.github.io/simIDM/reference/survPFSOS.md)
  : Survival Function of the Product PFS\*OS for Different Transition
  Models
- [`survOS()`](https://insightsengineering.github.io/simIDM/reference/survOS.md)
  : OS Survival Function for Different Transition Models
