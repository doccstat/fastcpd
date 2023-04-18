
<!-- README.md is generated from README.Rmd. Please edit that file -->

# fastcpd

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/fastcpd)](https://cran.r-project.org/package=fastcpd)
[![R-CMD-check](https://github.com/doccstat/fastcpd/workflows/R-CMD-check/badge.svg)](https://github.com/doccstat/fastcpd/actions)
[![Codecov test
coverage](https://codecov.io/gh/doccstat/fastcpd/branch/main/graph/badge.svg)](https://app.codecov.io/gh/doccstat/fastcpd?branch=main)
[![Last
Commit](https://img.shields.io/github/last-commit/doccstat/fastcpd)](https://github.com/doccstat/fastcpd)
<!-- badges: end -->

## Overview

The fastcpd (**FAST** \_\_C\_\_hange \_\_P\_\_oint \_\_D\_\_etection) is
a fast implmentation of change point detection methods in R. The
**fastcpd** package is designed to find change points in a fast manner.
It is easy to install and extensible to all kinds of change point
problems with a user specified cost function apart from the built-in
cost functions.

If you’d like to learn how to use the fastcpd effectively, please refer
to the following references:

- [Sequential Gradient Descent and Quasi-Newton’s Method for
  Change-Point
  Analysis](https://proceedings.mlr.press/v206/zhang23b.html)

## Installation

``` r
# Install from CRAN, (not yet available)
install.packages("fastcpd")
```

``` r
# Install the development version from GitHub
# install.packages("pak")
pak::pak("doccstat/fastcpd")
```

If you’re compiling from source, you can run
`pak::pkg_system_requirements("fastcpd")`, to see the complete set of
system packages needed on your machine.

## Usage

`library(fastcpd)` will load the core fastcpd package in full or in
part:

- [Rcpp](https://github.com/RcppCore/Rcpp), for C++ source code
  compilation.
- [RcppArmadillo](https://github.com/RcppCore/RcppArmadillo), for fast
  linear algebra.
- [fastglm](https://github.com/jaredhuling/fastglm), for fast
  generalized linear models.
- [DescTools](https://github.com/AndriSignorell/DescTools), for
  Winsorizing Poisson data.
- [glmnet](https://glmnet.stanford.edu/), for penalized regression.
- [ggplot2](https://github.com/tidyverse/ggplot2), for data
  visualization.

### Binomial

``` r
# TODO: examples in roxygen docs and multivariate custom cost functions, pkgdown references

# Linear regression

# Logistic regression
library(fastcpd)
set.seed(1)
x <- matrix(rnorm(1500, 0, 1), ncol = 5)
theta <- rbind(rnorm(5, 0, 1), rnorm(5, 2, 1))
y <- c(
  rbinom(125, 1, 1 / (1 + exp(-x[1:125, ] %*% theta[1, ]))),
  rbinom(175, 1, 1 / (1 + exp(-x[126:300, ] %*% theta[2, ])))
)
result <- fastcpd(
  formula = y ~ . - 1,
  data = data.frame(y = y, x = x),
  family = "binomial",
  cp_only = FALSE
)
plot(result)
summary(result)
#> Call:
#> fastcpd(formula = y ~ . - 1, data = data.frame(y = y, x = x),
#>     family = "binomial", cp_only = FALSE)
#>
#> Residuals:
#>       Min        1Q    Median        3Q       Max
#> -14.09576  -1.07218  -1.00000   1.07353  35.39472
#>
#> Change points:
#> 126

# Poisson regression
library(fastcpd)
set.seed(1)
x <- matrix(rnorm(1500, 0, 1), ncol = 5)
theta <- rbind(rnorm(5, 0, 1), rnorm(5, 2, 1))
y <- c(
  rpois(125, exp(x[1:125, ] %*% theta[1, ])),
  rpois(175, exp(x[126:300, ] %*% theta[2, ]))
)
result <- fastcpd(
  formula = y ~ . - 1,
  data = data.frame(y = y, x = x),
  family = "poisson",
  cp_only = FALSE
)
plot(result)
summary(result)

# Penalized linear regression

# Custom cost function: mean shift
library(fastcpd)
set.seed(1)
p <- 1
data <- rbind(
  mvtnorm::rmvnorm(300, mean = rep(0, p), sigma = diag(100, p)),
  mvtnorm::rmvnorm(400, mean = rep(50, p), sigma = diag(100, p)),
  mvtnorm::rmvnorm(300, mean = rep(2, p), sigma = diag(100, p))
)

segment_count_guess <- 10
block_size <- max(floor(sqrt(nrow(data)) / (segment_count_guess + 1)), 2)
block_count <- ceiling(nrow(data) / block_size)
data_all_vars <- rep(0, block_count)
for (block_index in seq_len(block_count)) {
  data_all_vars[block_index] <- var(data[((block_index - 1) * block_size + 1):(min(block_index * block_size, nrow(data))), , drop = FALSE])
}
data_all_var <- mean(data_all_vars)

mean_loss <- function(data) {
  (norm(data, type = "F") ^ 2 - colSums(data) ^ 2 / nrow(data)) / 2 / data_all_var + nrow(data) / 2 * (log(data_all_var) + log(2 * pi))
}

mean_loss_result <- fastcpd(
  formula = ~ . - 1,
  data = data.frame(data),
  beta = (p + 1) * log(nrow(data)) / 2,
  p = p,
  cost = mean_loss
)

summary(mean_loss_result)

# Custom cost function: variance change
library(fastcpd)
set.seed(1)
p <- 1
data <- rbind.data.frame(
  mvtnorm::rmvnorm(300, mean = rep(0, p), sigma = diag(1, p)),
  mvtnorm::rmvnorm(400, mean = rep(0, p), sigma = diag(50, p)),
  mvtnorm::rmvnorm(300, mean = rep(0, p), sigma = diag(2, p))
)
data_all_mu <- colMeans(data)

var_loss <- function(data) {
  demeaned_data_norm <- norm(sweep(data, 2, data_all_mu), type = "F")
  nrow(data) * (1 + log(2 * pi) + log(demeaned_data_norm ^ 2 / nrow(data))) / 2
}

var_loss_result <- fastcpd(
  formula = ~ . - 1,
  data = data,
  beta = (p + 1) * log(nrow(data)) / 2,
  p = p,
  cost = var_loss
)

summary(var_loss_result)

# Custom cost function: mean shift and variance change
library(fastcpd)
set.seed(1)
p <- 1
data <- rbind.data.frame(
  mvtnorm::rmvnorm(300, mean = rep(0, p), sigma = diag(1, p)),
  mvtnorm::rmvnorm(400, mean = rep(10, p), sigma = diag(1, p)),
  mvtnorm::rmvnorm(300, mean = rep(0, p), sigma = diag(50, p)),
  mvtnorm::rmvnorm(300, mean = rep(0, p), sigma = diag(1, p)),
  mvtnorm::rmvnorm(400, mean = rep(10, p), sigma = diag(1, p)),
  mvtnorm::rmvnorm(300, mean = rep(10, p), sigma = diag(50, p))
)

meanvar_loss <- function(data) {
  nrow(data) * (1 + log(2 * pi) + log((colSums(data ^ 2) - colSums(data) ^ 2 / nrow(data)) / nrow(data))) / 2
}

meanvar_loss_result <- fastcpd(
  formula = ~ . - 1,
  data = data,
  beta = (2 * p + 1) * log(nrow(data)) / 2,
  p = 2 * p,
  cost = meanvar_loss
)

summary(meanvar_loss_result)

# Custom cost function: Huber loss
library(fastcpd)
set.seed(1)
n <- 400 + 300 + 400
p <- 3
x <- mvtnorm::rmvnorm(n, mean = rep(0, p), sigma = diag(p))
theta <- rbind(
  mvtnorm::rmvnorm(1, mean = rep(0, p), sigma = diag(p)),
  mvtnorm::rmvnorm(1, mean = rep(3, p), sigma = diag(p)),
  mvtnorm::rmvnorm(1, mean = rep(5, p), sigma = diag(p))
)
theta <- theta[rep(seq_len(3), c(400, 300, 400)), ]
y_true <- rowSums(x * theta)
factor <- c(
  2 * stats::rbinom(400, size = 1, prob = 0.95) - 1,
  2 * stats::rbinom(300, size = 1, prob = 0.95) - 1,
  2 * stats::rbinom(400, size = 1, prob = 0.95) - 1
)
y <- factor * y_true + stats::rnorm(n)
data <- cbind.data.frame(y, x)

huber_threshold <- 1

huber_loss <- function(data, theta) {
  residual <- data[, 1] - data[, -1, drop = FALSE] %*% theta
  indicator <- abs(residual) <= huber_threshold
  sum(residual ^ 2 / 2 * indicator + huber_threshold * (abs(residual) - huber_threshold / 2) * (1 - indicator))
}

huber_loss_gradient <- function(data, theta) {
  residual <- c(data[nrow(data), 1] - data[nrow(data), -1] %*% theta)
  if (abs(residual) <= huber_threshold) {
    - residual * data[nrow(data), -1]
  } else {
    - huber_threshold * sign(residual) * data[nrow(data), -1]
  }
}

huber_loss_hessian <- function(data, theta) {
  residual <- c(data[nrow(data), 1] - data[nrow(data), -1] %*% theta)
  if (abs(residual) <= huber_threshold) {
    outer(data[nrow(data), -1], data[nrow(data), -1])
  } else {
    0.01 * diag(length(theta))
  }
}

huber_regression_result <- fastcpd(
  formula = y ~ . - 1,
  data = data,
  beta = (p + 1) * log(n) / 2,
  cost = huber_loss,
  cost_gradient = huber_loss_gradient,
  cost_hessian = huber_loss_hessian
)

summary(huber_regression_result)
```
