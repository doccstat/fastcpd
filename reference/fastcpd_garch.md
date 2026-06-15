# Find change points efficiently in GARCH(\\p\\, \\q\\) models

`fastcpd_garch()` and `fastcpd.garch()` are wrapper functions of
[`fastcpd()`](https://fastcpd.xingchi.li/reference/fastcpd.md) to find
change points in GARCH(\\p\\, \\q\\) models. The function is similar to
[`fastcpd()`](https://fastcpd.xingchi.li/reference/fastcpd.md) except
that the data is by default a one-column matrix or univariate vector and
thus a formula is not required here.

## Usage

``` r
fastcpd_garch(data, order = c(0, 0), ...)

fastcpd.garch(data, order = c(0, 0), ...)
```

## Arguments

- data:

  A numeric vector, a matrix, a data frame or a time series object.

- order:

  A positive integer vector of length two specifying the order of the
  GARCH model.

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
# \donttest{
set.seed(1)
n <- 1501
sigma_2 <- rep(1, n + 1)
x <- rep(0, n + 1)
for (i in seq_len(750)) {
  sigma_2[i + 1] <- 20 + 0.8 * x[i]^2 + 0.1 * sigma_2[i]
  x[i + 1] <- rnorm(1, 0, sqrt(sigma_2[i + 1]))
}
for (i in 751:n) {
  sigma_2[i + 1] <- 1 + 0.1 * x[i]^2 + 0.5 * sigma_2[i]
  x[i + 1] <- rnorm(1, 0, sqrt(sigma_2[i + 1]))
}
result <- suppressWarnings(
  fastcpd.garch(x[-1], c(1, 1), include.mean = FALSE)
)
summary(result)
#> 
#> Call:
#> fastcpd.garch(data = x[-1], order = c(1, 1), include.mean = FALSE)
#> 
#> Change points:
#> 756 
#> 
#> Cost values:
#> 1959.888 742.4338 
#> 
#> Parameters:
#>    segment 1 segment 2
#> 1 19.0363193 1.5162583
#> 2  0.6394788 0.1329554
#> 3  0.2164783 0.2969624
plot(result)

# }
set.seed(1)
n <- 200
sigma_2 <- rep(1, n + 1)
x <- rep(0, n + 1)
for (i in seq_len(100)) {
  sigma_2[i + 1] <- 10 + 0.5 * x[i]^2 + 0.3 * sigma_2[i]
  x[i + 1] <- rnorm(1, 0, sqrt(sigma_2[i + 1]))
}
for (i in 101:n) {
  sigma_2[i + 1] <- 0.2 + 0.05 * x[i]^2 + 0.1 * sigma_2[i]
  x[i + 1] <- rnorm(1, 0, sqrt(sigma_2[i + 1]))
}
result <- suppressWarnings(
  fastcpd.garch(x[-1], c(1, 1), include.mean = FALSE)
)
summary(result)
#> 
#> Call:
#> fastcpd.garch(data = x[-1], order = c(1, 1), include.mean = FALSE)
#> 
#> Change points:
#> 98 
#> 
#> Cost values:
#> 202.3242 -16.34071 
#> 
#> Parameters:
#>   segment 1  segment 2
#> 1 5.6718832 0.19442986
#> 2 0.1830536 0.08163262
#> 3 0.5651263 0.04133033
plot(result)

set.seed(1)
n <- 300
sigma_2 <- rep(1, n + 1)
x <- rep(0, n + 1)
for (i in seq_len(100)) {
  sigma_2[i + 1] <- 10 + 0.5 * x[i]^2 + 0.3 * sigma_2[i]
  x[i + 1] <- rnorm(1, 0, sqrt(sigma_2[i + 1]))
}
for (i in 101:n) {
  sigma_2[i + 1] <- 0.2 + 0.05 * x[i]^2 + 0.1 * sigma_2[i]
  x[i + 1] <- rnorm(1, 0, sqrt(sigma_2[i + 1]))
}
result <- suppressWarnings(
  fastcpd.garch(x[-1], c(1, 1), include.mean = FALSE)
)
summary(result)
#> 
#> Call:
#> fastcpd.garch(data = x[-1], order = c(1, 1), include.mean = FALSE)
#> 
#> Change points:
#> 98 
#> 
#> Cost values:
#> 202.3242 -34.68556 
#> 
#> Parameters:
#>   segment 1  segment 2
#> 1 5.6718832 0.21471834
#> 2 0.1830536 0.09268943
#> 3 0.5651263 0.00000000
plot(result)
```
