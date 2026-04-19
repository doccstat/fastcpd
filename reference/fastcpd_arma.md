# Find change points efficiently in ARMA(\\p\\, \\q\\) models

`fastcpd_arma()` and `fastcpd.arma()` are wrapper functions of
[`fastcpd()`](https://fastcpd.xingchi.li/reference/fastcpd.md) to find
change points in ARMA(\\p\\, \\q\\) models. The function is similar to
[`fastcpd()`](https://fastcpd.xingchi.li/reference/fastcpd.md) except
that the data is by default a one-column matrix or univariate vector and
thus a formula is not required here.

## Usage

``` r
fastcpd_arma(data, order = c(0, 0), ...)

fastcpd.arma(data, order = c(0, 0), ...)
```

## Arguments

- data:

  A numeric vector, a matrix, a data frame or a time series object.

- order:

  A vector of length two specifying the order of the ARMA model.

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
n <- 200
w <- rnorm(n + 3, 0, 3)
x <- rep(0, n + 3)
for (i in 1:150) {
  x[i + 3] <- 0.1 * x[i + 2] - 0.3 * x[i + 1] + 0.1 * x[i] +
    0.1 * w[i + 2] + 0.5 * w[i + 1] + w[i + 3]
}
for (i in 151:n) {
  x[i + 3] <- 0.3 * x[i + 2] + 0.1 * x[i + 1] - 0.3 * x[i] -
    0.6 * w[i + 2] - 0.1 * w[i + 1] + w[i + 3]
}
result <- suppressWarnings(
  fastcpd.arma(
    data = x[3 + seq_len(n)],
    order = c(3, 2),
    segment_count = 3,
    lower = c(rep(-1, 3 + 2), 1e-10),
    upper = c(rep(1, 3 + 2), Inf),
    line_search = c(1, 0.1, 1e-2),
    beta = "BIC",
    cost_adjustment = "BIC",
    trim = 0.025
  )
)
summary(result)
#> 
#> Call:
#> fastcpd.arma(data = x[3 + seq_len(n)], order = c(3, 2), segment_count = 3, 
#>     lower = c(rep(-1, 3 + 2), 1e-10), upper = c(rep(1, 3 + 2), 
#>         Inf), line_search = c(1, 0.1, 0.01), beta = "BIC", cost_adjustment = "BIC", 
#>     trim = 0.025)
#> 
#> Change points:
#> 139 
#> 
#> Cost values:
#> 330.6779 157.0954 
#> 
#> Parameters:
#>    segment 1  segment 2
#> 1  0.7108329  1.2165407
#> 2 -0.1568921 -0.7959726
#> 3 -0.1388433 -0.1341928
#> 4 -0.5309670 -1.4286676
#> 5  0.3217713  0.9999122
#> 6  6.8097404  9.4577537
plot(result)

# }
```
