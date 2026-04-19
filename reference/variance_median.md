# Variance estimation for median change models

Implement Rice estimator.

## Usage

``` r
variance_median(data)

variance.median(data)
```

## Arguments

- data:

  A vector of data points.

## Value

A numeric value representing the variance.

## Examples

``` r
(sigma2 <- variance.median(well_log))
#> [1] 5803645
```
