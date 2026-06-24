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
n <- 401
sigma_2 <- rep(1, n + 1)
x <- rep(0, n + 1)
for (i in seq_len(200)) {
  sigma_2[i + 1] <- 20 + 0.8 * x[i]^2 + 0.1 * sigma_2[i]
  x[i + 1] <- rnorm(1, 0, sqrt(sigma_2[i + 1]))
}
for (i in 201:n) {
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
#> 205 
#> 
#> Cost values:
#> 491.9506 188.9431 
#> 
#> Parameters:
#>    segment 1 segment 2
#> 1 11.5363695 1.9983858
#> 2  0.4695383 0.0128258
#> 3  0.3428052 0.1368199
plot(result)

# }
# \donttest{
set.seed(1)
n <- 120
sigma_2 <- rep(1, n + 1)
x <- rep(0, n + 1)
for (i in seq_len(60)) {
  sigma_2[i + 1] <- 10 + 0.5 * x[i]^2 + 0.3 * sigma_2[i]
  x[i + 1] <- rnorm(1, 0, sqrt(sigma_2[i + 1]))
}
for (i in 61:n) {
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
#> 60 
#> 
#> Cost values:
#> 120.2821 -15.87097 
#> 
#> Parameters:
#>    segment 1 segment 2
#> 1 13.8346334  0.174438
#> 2  0.2278886  0.000000
#> 3  0.0200848  0.000000
plot(result)

# }
# \donttest{
set.seed(1)
n <- 150
sigma_2 <- rep(1, n + 1)
x <- rep(0, n + 1)
for (i in seq_len(75)) {
  sigma_2[i + 1] <- 10 + 0.5 * x[i]^2 + 0.3 * sigma_2[i]
  x[i + 1] <- rnorm(1, 0, sqrt(sigma_2[i + 1]))
}
for (i in 76:n) {
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
#> 74 
#> 
#> Cost values:
#> 154.8803 -19.62733 
#> 
#> Parameters:
#>   segment 1   segment 2
#> 1 2.8412668 0.180061756
#> 2 0.1358175 0.005695202
#> 3 0.7489304 0.000000000
plot(result)

# }
```
