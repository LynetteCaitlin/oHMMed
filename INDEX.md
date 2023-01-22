
<!-- INDEX.md is generated from INDEX.Rmd. Please edit that file -->

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
devtools::install_github("majkamichal/oHMMed")
```
