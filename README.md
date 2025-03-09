
<!-- README.md is generated from README.Rmd. Please edit that file -->

# fastcpd: Fast Change Point Detection in R <a href="https://fastcpd.xingchi.li"><img src="man/figures/logo.svg" align="right" height="138" /></a>

[![Codecov test
coverage](https://codecov.io/gh/doccstat/fastcpd/branch/main/graph/badge.svg)](https://app.codecov.io/gh/doccstat/fastcpd?branch=main)
[![CodeFactor](https://www.codefactor.io/repository/github/doccstat/fastcpd/badge)](https://www.codefactor.io/repository/github/doccstat/fastcpd)
[![CRAN
status](https://www.r-pkg.org/badges/version-last-release/fastcpd)](https://cran.r-project.org/package=fastcpd)
[![doi](https://img.shields.io/badge/doi-10.48550/arXiv.2404.05933-green.svg)](https://doi.org/10.48550/arXiv.2404.05933)
[![Downloads](https://cranlogs.r-pkg.org/badges/grand-total/fastcpd)](https://cran.r-project.org/package=fastcpd)
[![Last
Commit](https://img.shields.io/github/last-commit/doccstat/fastcpd)](https://github.com/doccstat/fastcpd)
[![R-CMD-check](https://github.com/doccstat/fastcpd/workflows/R-CMD-check/badge.svg)](https://github.com/doccstat/fastcpd/actions)
[![r-universe](https://doccstat.r-universe.dev/badges/fastcpd)](https://doccstat.r-universe.dev)

# Branch `profiling` is for profiling purposes only.

``` r
library(fastcpd)
set.seed(1)
p <- 3
data <- rbind(
  mvtnorm::rmvnorm(3e+6, mean = rep(0, p), sigma = diag(100, p)),
  mvtnorm::rmvnorm(4e+6, mean = rep(50, p), sigma = diag(100, p)),
  mvtnorm::rmvnorm(3e+6, mean = rep(2, p), sigma = diag(100, p))
)

fastcpd:::start_profiler("./profile.out")
fastcpd.mean(data, r.progress = FALSE, cp_only = TRUE)
fastcpd:::stop_profiler()
```

``` shell
pprof --pdf /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/fastcpd/libs/fastcpd.so ./profile.out > profile.pdf
```

## Overview

The fastcpd (**fast** **c**hange **p**oint **d**etection) is a fast
implmentation of change point detection methods in R. The **fastcpd**
package is designed to find change points in a fast manner. It is easy
to install and extensible to all kinds of change point problems with a
user specified cost function apart from the built-in cost functions.

To learn more behind the algorithms:

- [fastcpd: Fast Change Point Detection in
  R](https://doi.org/10.48550/arXiv.2404.05933)
- [Sequential Gradient Descent and Quasi-Newton’s Method for
  Change-Point
  Analysis](https://proceedings.mlr.press/v206/zhang23b.html)

## Installation

``` r
# Install from r-universe with CRAN version as a fallback
install.packages(
  "fastcpd",
  repos = c("https://doccstat.r-universe.dev", "https://cloud.r-project.org")
)
```

<details close>
<summary>
Development version
</summary>

``` r
## install.packages("pak")
pak::pak("doccstat/fastcpd")

## install.packages("devtools")
devtools::install_github("doccstat/fastcpd")
```

</details>
<details close>
<summary>
With mamba or conda (available soon)
</summary>

``` bash
# conda-forge is a fork from CRAN and may not be up-to-date

# Use mamba
mamba install r-fastcpd
# Use conda
conda install -c conda-forge r-fastcpd
```

</details>

### FAQ

<details close>
<summary>
Package dependencies
</summary>

`fastcpd` depends on the following packages:

- [Rcpp](https://github.com/RcppCore/Rcpp), for C++ source code
  compilation.
- [RcppArmadillo](https://github.com/RcppCore/RcppArmadillo), for fast
  linear algebra.
- [fastglm](https://github.com/jaredhuling/fastglm), for fast
  generalized linear models.
- [glmnet](https://glmnet.stanford.edu/), for penalized regression.
- [ggplot2](https://github.com/tidyverse/ggplot2), for data
  visualization.

If you’re compiling from source, you can run the following command to
see the complete set of system packages needed on your machine.

``` r
pak::pkg_sysreqs("doccstat/fastcpd")
#> ── Install scripts ───────────────────────────────────────────── Ubuntu 20.04
#> apt-get -y update
#> apt-get -y install libcurl4-openssl-dev libssl-dev zlib1g-dev make
#>
#> ── Packages and their system dependencies ───────────────────────────────────
#> curl       – libcurl4-openssl-dev, libssl-dev
#> data.table – zlib1g-dev
#> fs         – make
#> openssl    – libssl-dev
```

</details>
<details close>
<summary>
Should I install suggested packages?
</summary>

The suggested packages are not required for the main functionality of
the package. They are only required for the vignettes. If you want to
learn more about the package comparison and other vignettes, you could
either check out vignettes on
[CRAN](https://CRAN.R-project.org/package=fastcpd) or [pkgdown generated
documentation](https://fastcpd.xingchi.li/articles/).

</details>
<details close>
<summary>
I countered problems related to gfortran on Mac OSX or Linux!
</summary>

The package should be able to install on Mac and any Linux distribution
without any problems if all the dependencies are installed. However, if
you encountered problems related to gfortran, it might be because
`RcppArmadillo` is not installed previously. Try [Mac OSX stackoverflow
solution](https://stackoverflow.com/a/72997915) or [Linux stackover
solution](https://stackoverflow.com/a/15540919) if you have trouble
installing `RcppArmadillo`.

</details>

## Cheatsheet

[![fastcpd
cheatsheet](man/figures/cheatsheets.png)](https://github.com/doccstat/fastcpd/blob/main/man/figures/cheatsheets.pdf)

### R Shiny App

Available soon:
[rshiny.fastcpd.xingchi.li](https://rshiny.fastcpd.xingchi.li)

## Usage

``` r
set.seed(1)
n <- 1000
x <- rep(0, n + 3)
for (i in 1:600) {
  x[i + 3] <- 0.6 * x[i + 2] - 0.2 * x[i + 1] + 0.1 * x[i] + rnorm(1, 0, 3)
}
for (i in 601:1000) {
  x[i + 3] <- 0.3 * x[i + 2] + 0.4 * x[i + 1] + 0.2 * x[i] + rnorm(1, 0, 3)
}
result <- fastcpd::fastcpd.ar(x[3 + seq_len(n)], 3, r.progress = FALSE)
summary(result)
#> 
#> Call:
#> fastcpd::fastcpd.ar(data = x[3 + seq_len(n)], order = 3, r.progress = FALSE)
#> 
#> Change points:
#> 614 
#> 
#> Cost values:
#> 2754.116 2038.945 
#> 
#> Parameters:
#>     segment 1 segment 2
#> 1  0.57120256 0.2371809
#> 2 -0.20985108 0.4031244
#> 3  0.08221978 0.2290323
plot(result)
```

![](man/figures/README-ar3-1.png)<!-- -->

<details open>
<summary>
Tip
</summary>

It is hard to demonstrate all the features of `fastcpd` in a single
example due to the flexibility of the package. For more examples, please
refer to the [function
reference](https://fastcpd.xingchi.li/reference/index.html).

</details>
<details open>
<summary>
Note
</summary>

`r.progress = FALSE` is used to suppress the progress bar. Users are
expected to see the progress bar when running the code by default.

</details>

## Examples

### Main function

- [`fastcpd`](https://fastcpd.xingchi.li/reference/fastcpd.html)

### Wrapper functions

#### Time series

- AR(p):
  [`fastcpd_ar`](https://fastcpd.xingchi.li/reference/fastcpd_ar.html)
- ARIMA(p, d, q):
  [`fastcpd_arima`](https://fastcpd.xingchi.li/reference/fastcpd_arima.html)
- ARMA(p, q):
  [`fastcpd_arma`](https://fastcpd.xingchi.li/reference/fastcpd_arma.html)
- GARCH(p, q):
  [`fastcpd_garch`](https://fastcpd.xingchi.li/reference/fastcpd_garch.html)
- VAR(p):
  [`fastcpd_var`](https://fastcpd.xingchi.li/reference/fastcpd_var.html)
- General time series:
  [`fastcpd_ts`](https://fastcpd.xingchi.li/reference/fastcpd_ts.html)

#### Unlabeled data

- Mean change:
  [`fastcpd_mean`](https://fastcpd.xingchi.li/reference/fastcpd_mean.html)
- Variance change:
  [`fastcpd_variance`](https://fastcpd.xingchi.li/reference/fastcpd_variance.html)
- Mean and/or variance change:
  [`fastcpd_meanvariance`](https://fastcpd.xingchi.li/reference/fastcpd_meanvariance.html)

#### Regression data

- Logistic regression:
  [`fastcpd_binomial`](https://fastcpd.xingchi.li/reference/fastcpd_binomial.html)
- Penalized linear regression:
  [`fastcpd_lasso`](https://fastcpd.xingchi.li/reference/fastcpd_lasso.html)
- Linear regression:
  [`fastcpd_lm`](https://fastcpd.xingchi.li/reference/fastcpd_lm.html)
- Poisson regression:
  [`fastcpd_poisson`](https://fastcpd.xingchi.li/reference/fastcpd_poisson.html)

### Utility functions

#### Variance estimation

- Variance estimation in ARMA models:
  [`variance_arma`](https://fastcpd.xingchi.li/reference/variance_arma.html)
- Variance estimation in linear models:
  [`variance_lm`](https://fastcpd.xingchi.li/reference/variance_lm.html)
- Variance estimation in mean change models:
  [`variance_mean`](https://fastcpd.xingchi.li/reference/variance_mean.html)
- Variance estimation in median change models:
  [`variance_median`](https://fastcpd.xingchi.li/reference/variance_median.html)

#### Class methods

- Plot: [`plot`](https://fastcpd.xingchi.li/reference/plot.html)
- Print: [`print`](https://fastcpd.xingchi.li/reference/print.html)
- Show: [`show`](https://fastcpd.xingchi.li/reference/show.html)
- Summary:
  [`summary`](https://fastcpd.xingchi.li/reference/summary.html)

### Data

- Bitcoin Market Price (USD):
  [`bitcoin`](https://fastcpd.xingchi.li/reference/bitcoin.html)
- Occupancy Detection Data Set:
  [`occupancy`](https://fastcpd.xingchi.li/reference/occupancy.html)
- Transcription Profiling of 57 Human Bladder Carcinoma Samples:
  [`transcriptome`](https://fastcpd.xingchi.li/reference/transcriptome.html)
- UK Seatbelts Data:
  [`uk_seatbelts`](https://fastcpd.xingchi.li/reference/uk_seatbelts.html)
- Well-log Dataset from Numerical Bayesian Methods Applied to Signal
  Processing:
  [`well_log`](https://fastcpd.xingchi.li/reference/well_log.html)

### Main class

- [`fastcpd`](https://fastcpd.xingchi.li/reference/fastcpd-class.html)

## Make contributions

<details close>
<summary>
We welcome contributions from everyone. Please follow the instructions
below to make contributions.
</summary>

1.  Fork the repo.

2.  Create a new branch from `main` branch.

3.  Make changes and commit them.

    1.  Please follow the [Google’s R style
        guide](https://google.github.io/styleguide/Rguide.html) for
        naming variables and functions.
    2.  If you are adding a new family of models with new cost functions
        with corresponding gradient and Hessian, please add them to
        `src/fastcpd_class_cost.cc` with proper example and tests in
        `vignettes/gallery.Rmd` and `tests/testthat/test-gallery.R`.
    3.  Add the family name to `src/fastcpd_constants.h`.
    4.  \[Recommended\] Add a new wrapper function in
        `R/fastcpd_wrappers.R` for the new family of models and move the
        examples to the new wrapper function as roxygen examples.
    5.  Add the new wrapper function to the corresponding section in
        `_pkgdown.yml`.

4.  Push the changes to your fork.

5.  Create a pull request.

6.  Make sure the pull request does not create new warnings or errors in
    `devtools::check()`.

</details>

## Contact us

<details close>
<summary>
Encountered a bug or unintended behavior?
</summary>

1.  File a ticket at [GitHub
    Issues](https://github.com/doccstat/fastcpd/issues).
2.  Contact the authors specified in
    [DESCRIPTION](https://github.com/doccstat/fastcpd/blob/main/DESCRIPTION#L5-L10).

</details>

## Stargazers over time

[![Stargazers over
time](https://starchart.cc/doccstat/fastcpd.svg)](https://starchart.cc/doccstat/fastcpd)

## Codecov Icicle

[![Codecov
Icicle](https://codecov.io/gh/doccstat/fastcpd/graphs/icicle.svg)](https://app.codecov.io/gh/doccstat/fastcpd?branch=main)
