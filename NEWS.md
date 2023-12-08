# simIDM 0.1.0

### New Features

-   `estimateParams` estimates parameters of exponential and Weibull transition hazards using maximum likelihood for a given data set after using `prepareData`.
-   `corTrans` calculates the correlation between PFS and OS for a given set of transition hazards, and `corPFSOS` estimates this correlation for a given data set, optionally with a bootstrap based confidence interval.

### Bug Fixes

-   `ExpSurvOS` now returns 0 instead of NaN for large values of t.
-   `WeibSurvOS` now does not return an error for large values of t.
-   `PWCSurvOS` now does not return an error for large values of t. It also no longer returns values larger than 1. It is significantly faster, based on a closed form calculation instead of numerical integration.
-   `getSimulatedData` now also works when there are no transitions from progression to death, similarly for `getOneClinicalTrial` (which now warns if there are no such transitions at all).

### Miscellaneous

-   `PwcOSInt`, `integrateVector`, `WeibOSInteg` are no longer exported, and only used internally.
-   Renamed piecewise constant hazards function to `getPWCHazard` (previously `getPCWHazard`).

# simIDM 0.0.5

-   First CRAN version of the package.
-   The package simulates illness-death models with constant, Weibull or piecewise constant transition hazards.

### New Features

-   Exponentially, Weibull and piecewise exponentially distributed survival times.
-   Random censoring and event-driven censoring after a pre-specified number of PFS or OS events.
-   Arbitrary number of treatment arms and flexible randomization ratio.
-   Staggered study entry.
-   Derivation of PFS and OS survival functions from transition hazards.
