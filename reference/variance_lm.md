# Variance estimation for linear models with change points

Estimate the variance for each block and then take the average.

## Usage

``` r
variance_lm(data, d = 1, block_size = ncol(data) - d + 1, outlier_iqr = Inf)

variance.lm(data, d = 1, block_size = ncol(data) - d + 1, outlier_iqr = Inf)
```

## Arguments

- data:

  A matrix or a data frame with the response variable as the first
  column.

- d:

  The dimension of the response variable.

- block_size:

  The size of the blocks to use for variance estimation.

- outlier_iqr:

  The number of interquartile ranges to use as a threshold for outlier
  detection.

## Value

A numeric value representing the variance.

## Examples

``` r
if (requireNamespace("mvtnorm", quietly = TRUE)) {
  set.seed(1)
  n <- 300
  p <- 4
  x <- mvtnorm::rmvnorm(n, rep(0, p), diag(p))
  theta <- rbind(c(1, 3.2, -1, 0), c(-1, -0.5, 2.5, -2), c(0.8, 0, 1, 2))
  y <- c(
    x[1:100, ] %*% theta[1, ] + rnorm(100, 0, 3),
    x[101:200, ] %*% theta[2, ] + rnorm(100, 0, 3),
    x[201:n, ] %*% theta[3, ] + rnorm(100, 0, 3)
  )
  (sigma2 <- variance.lm(cbind(y, x)))

  set.seed(1)
  n <- 300
  p <- 4
  d <- 2
  x <- mvtnorm::rmvnorm(n, rep(0, p), diag(p))
  theta <- cbind(c(1, 3.2, -1, 0), c(-1, -0.5, 2.5, -2), c(0.8, 0, 1, 2))
  theta <- cbind(theta, theta)
  y <- rbind(
    x[1:100, ] %*% theta[, 1:2] +
      mvtnorm::rmvnorm(100, rep(0, d), diag(3, d)),
    x[101:200, ] %*% theta[, 3:4] +
      mvtnorm::rmvnorm(100, rep(0, d), diag(3, d)),
    x[201:n, ] %*% theta[, 5:6] +
      mvtnorm::rmvnorm(100, rep(0, d), diag(3, d))
  )
  (sigma <- variance.lm(cbind(y, x), d = d))
}
#>            [,1]       [,2]
#> [1,] 3.71946553 0.01962422
#> [2,] 0.03459023 3.59825273
```
