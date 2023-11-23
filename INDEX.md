
<!-- INDEX.md is generated from INDEX.Rmd. Please edit that file -->

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/oHMMed)](https://cran.r-project.org/package=oHMMed)
[![](http://cranlogs.r-pkg.org/badges/oHMMed)](http://cran.rstudio.com/web/packages/oHMMed/index.html)
[![](https://cranlogs.r-pkg.org/badges/grand-total/oHMMed)](http://cran.rstudio.com/web/packages/oHMMed/index.html)

<br>

# oHMMed

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
parametrised so that they can be ordered by increasing mean. In fact,
transitions between states in the hidden sequence can only occur between
states that emit densities that are neighbours in the ordering by mean.
Given the observed data sequence, our models assign each part of the
hidden sequence to a state, and infer the transition rates between them
as well as the parameters of the state-specific emission distributions.
This is the general framework of `oHMMed` (ordered Hidden Markov Model
with emission densities), and it can be applied to any system that
fulfils these assumptions.

------------------------------------------------------------------------

## 2. Installation

You can install the latest stable
[cran](https://cran.r-project.org/web/packages/oHMMed/index.html)
version using (recommended):

``` r
install.packages("oHMMed")
```

In order to install the latest stable development version from GitHub
you can use:

``` r
# install.packages("devtools")
devtools::install_github("LynetteCaitlin/oHMMed@*release")
```

The development version can be installed using:

``` r
# install.packages("devtools")
devtools::install_github("LynetteCaitlin/oHMMed")
```

------------------------------------------------------------------------

## 3. Usage Recommendations

Usage recommendations can be found in [this
file](https://github.com/LynetteCaitlin/oHMMed/blob/main/UsageRecommendations.pdf).

<br>
