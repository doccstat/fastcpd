# Variance estimation for mean change models

Implement Rice estimator for variance in mean change models.

## Usage

``` r
variance_mean(data)

variance.mean(data)
```

## Arguments

- data:

  A matrix or a data frame with data points as each row.

## Value

A matrix representing the variance-covariance matrix or a numeric value
representing the variance.

## Examples

``` r
if (requireNamespace("mvtnorm", quietly = TRUE)) {
  set.seed(1)
  p <- 3
  data <- rbind(
    mvtnorm::rmvnorm(300, mean = rep(0, p), sigma = diag(100, p)),
    mvtnorm::rmvnorm(400, mean = rep(50, p), sigma = diag(100, p)),
    mvtnorm::rmvnorm(300, mean = rep(2, p), sigma = diag(100, p))
  )
  (sigma <- variance.mean(data))
}
#>             [,1]        [,2]       [,3]
#> [1,] 113.5442077  -0.1433207  -2.973472
#> [2,]  -0.1433207 106.1627730   8.430343
#> [3,]  -2.9734719   8.4303433 112.602989
```
