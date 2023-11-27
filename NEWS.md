# simIDM 0.0.5.9020

### New Features

-   `prepareData` allows formatting of trial data for log-likelihood computation.
-   `estimateParams` estimates parameters of exponential and Weibull transition hazards for a given data set.

### Bug Fixes

-   `ExpSurvOS` now returns 0 instead of NaN for large values of t.
-   `WeibSurvOS` now does not return an error for large values of t.
-   `PWCSurvOS` now does not return an error for large values of t.
-   `getSimulatedData` now also works when there are no transitions from progression to death, similarly for `getOneClinicalTrial` (which now warns if there are no such transitions at all).

# simIDM 0.0.5

-   First CRAN version of the package.
-   The package simulates illness-death models with constant, Weibull or piecewise constant transition hazards.

### New Features

-   Exponentially, Weibull and piecewise exponentially distributed survival times.
-   Random censoring and event-driven censoring after a pre-specified number of PFS or OS events.
-   Arbitrary number of treatment arms and flexible randomization ratio.
-   Staggered study entry.
-   Derivation of PFS and OS survival functions from transition hazards.
