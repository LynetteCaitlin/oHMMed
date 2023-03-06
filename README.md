
<!-- README.md is generated from README.Rmd. Please edit that file -->

Extended documentation can be found on the website:
<https://lynettecaitlin.github.io/oHMMed/>
(Website is still work in progress!)

## 1. Overview

The `oHMMed` package contains an implementation of Hidden Markov Models
with ordered hidden states and emission densities. More precisely: We
assume a sequence of un-observable (’hidden’) variables with discrete
categories known as states. Moving along this sequence, the probability
of a specific state occurring at any step depends only on the state
preceding it. Each hidden variable along the sequence produces/’emits’
an observable data point. We assume that these emitted data points are
distributed according to some continuous distribution (currently: normal
or gamma-poisson compound) with state-specific parameters. Further, we
assume that the continuous emission distributions per state are
parameterised so that they can be ordered by increasing mean. In fact,
transitions between states in the hidden sequence can only occur between
states that emit densities that are neighbours in the ordering by mean.
Given the observed data sequence, our models assign each part of the
hidden sequence to a state, and infer the transition rates between them
as well as the parameters of the state-specific emission distributions.
This is the general framework of `oHMMed` (ordered Hidden Markov Model
with emission densities), and it can be applied to any system that
fulfills these assumptions.

The algorithms are from the following paper: (...tba...)

## 2. Installation

<!-- Just like many other `R` packages, `oHMMed` can be installed from the `CRAN` repository by simply executing in the console the following line: -->
<!-- ```{r, eval = FALSE} -->
<!-- # install.packages("oHMMed") -->
<!-- # Or the the development version from GitHub: -->
<!-- devtools::install_github("majkamichal/oHMMed") -->
<!-- ``` -->

The `oHMMed` package is currently in development and it can be installed
from `github` repository by simply executing in the console the
following line:

``` r
devtools::install_github("LynetteCaitlin/oHMMed")
```

## 3. Usage
Please read the following usage recommendations: <https://github.com/LynetteCaitlin/oHMMed/blob/main/UsageRecommendations.pdf>

The accompanying R scripts are, for the case of normal emission densities...:

<https://github.com/LynetteCaitlin/oHMMed/blob/main/simulation_scripts/MCMC_normal_sims.R>

... and for gamma-poisson emission densities:

<https://github.com/LynetteCaitlin/oHMMed/blob/main/simulation_scripts/MCMC_gamma-poisson_sims.R>
 
Template code viewing the results of two different sequence annotations side-by-side and assessing associations between the two observed sequences can be found at:

<https://github.com/LynetteCaitlin/oHMMed/blob/main/simulation_scripts/StandardOutputPlots.R>

## 4. Package Maintainers

Michal Majka (R technicalities and code optimisation)

Lynette Caitlin Mikula (Algorithm development, usage, and application)

