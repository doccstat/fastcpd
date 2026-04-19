# Variance estimation for ARMA model with change points

Estimate the variance for each block and then take the average.

## Usage

``` r
variance_arma(data, p, q, max_order = p * q)

variance.arma(data, p, q, max_order = p * q)
```

## Arguments

- data:

  A one-column matrix or a vector.

- p:

  The order of the autoregressive part.

- q:

  The order of the moving average part.

- max_order:

  The maximum order of the AR model to consider.

## Value

A numeric value representing the variance.

## Examples

``` r
set.seed(1)
n <- 300
w <- rnorm(n + 3, 0, 10)
x <- rep(0, n + 3)
for (i in 1:200) {
  x[i + 3] <- 0.1 * x[i + 2] - 0.3 * x[i + 1] + 0.1 * x[i] +
    0.1 * w[i + 2] + 0.5 * w[i + 1] + w[i + 3]
}
for (i in 201:n) {
  x[i + 3] <- 0.3 * x[i + 2] + 0.1 * x[i + 1] - 0.3 * x[i] -
    0.6 * w[i + 2] - 0.1 * w[i + 1] + w[i + 3]
}
(result <- variance.arma(x[-seq_len(3)], p = 3, q = 2))
#> $table
#>         sigma2      AIC      BIC
#> AR(1) 104.8647 4.659337 4.671683
#> AR(2) 100.2594 4.621094 4.645786
#> AR(3) 115.3102 4.767626 4.804664
#> AR(4) 105.4397 4.684806 4.734190
#> AR(5) 115.4258 4.781962 4.843691
#> AR(6) 117.0668 4.802745 4.876821
#> 
#> $sigma2_aic
#> [1] 100.2594
#> 
#> $sigma2_bic
#> [1] 100.2594
#> 
```
