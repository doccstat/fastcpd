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
#> 759 
#> 
#> Cost values:
#> 2034.494 736.2013 
#> 
#> Parameters:
#>      segment 1 segment 2
#> 1 7.982285e+01 1.3761900
#> 2 1.965442e-01 0.1254212
#> 3 7.275673e-13 0.3565612
plot(result)

# }
```
