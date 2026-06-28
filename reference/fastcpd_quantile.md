# Find change points efficiently in quantile regression models

`fastcpd_quantile()` and `fastcpd.quantile()` are wrapper functions of
[`fastcpd()`](https://fastcpd.xingchi.li/reference/fastcpd.md) to detect
change points in quantile regression models using the pinball (check
function) loss \\\rho\_\tau(u) = u(\tau - \mathbf{1}\_{u \< 0})\\. The
function detects changes in the conditional \\\tau\\-quantile of the
response given the covariates. The segment cost is minimised via
iteratively reweighted least squares (IRLS).

## Usage

``` r
fastcpd_quantile(data, order = 0.5, ...)

fastcpd.quantile(data, order = 0.5, ...)
```

## Arguments

- data:

  A matrix or a data frame with the response variable as the first
  column and covariates in the remaining columns.

- order:

  Quantile level \\\tau\\, a numeric value in (0, 1). The default is
  0.5, corresponding to median regression.

- ...:

  Other arguments passed to
  [`fastcpd()`](https://fastcpd.xingchi.li/reference/fastcpd.md), for
  example, `segment_count`.

## Value

A [fastcpd](https://fastcpd.xingchi.li/reference/fastcpd-class.md)
object.

## See also

[`fastcpd()`](https://fastcpd.xingchi.li/reference/fastcpd.md)

## Examples

``` r
set.seed(1)
n <- 200
p <- 2
x <- matrix(rnorm(n * p), n, p)
theta_1 <- c(1, -1)
theta_2 <- c(-1, 1)
y <- c(
  x[1:100, ] %*% theta_1 + rnorm(100),
  x[101:200, ] %*% theta_2 + rnorm(100)
)
result <- fastcpd_quantile(cbind(y, x), order = 0.5)
summary(result)
#> 
#> Call:
#> fastcpd_quantile(data = cbind(y, x), order = 0.5)
#> 
#> Change points:
#> 100 
#> 
#> Cost values:
#> 50.56047 43.33312 
#> 
#> Parameters:
#>    segment 1  segment 2
#> 1  0.9882078 -1.1610592
#> 2 -1.0449301  0.9042203
set.seed(42)
n <- 400
x <- rep(1, n)
eps <- rnorm(n)
eps[c(50, 51, 52, 300, 301, 302)] <- -15
y <- c(0 + eps[1:200], 3 + eps[201:400])
dat <- cbind(y, x)
result_lm <- fastcpd_lm(dat)
summary(result_lm)
#> 
#> Call:
#> fastcpd_lm(data = dat)
#> 
#> Change points:
#> 49 50 51 52 201 299 300 301 302 
#> 
#> Cost values:
#> 34.18657 0 0 0 64.69095 50.13901 0 0 0 40.02713 
#> 
#> Parameters:
#>     segment 1 segment 2 segment 3 segment 4   segment 5 segment 6 segment 7
#> 1 -0.04978035       -15       -15       -15 -0.01511621  3.026426       -12
#>   segment 8 segment 9 segment 10
#> 1       -12       -12   3.025898
result_qr <- fastcpd_quantile(dat, order = 0.5)
summary(result_qr)
#> 
#> Call:
#> fastcpd_quantile(data = dat, order = 0.5)
#> 
#> Change points:
#> 201 
#> 
#> Cost values:
#> 101.4315 97.94862 
#> 
#> Parameters:
#>     segment 1 segment 2
#> 1 -0.02944488  2.946003
```
