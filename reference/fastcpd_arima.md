# Find change points efficiently in ARIMA(\\p\\, \\d\\, \\q\\) models

`fastcpd_arima()` and `fastcpd.arima()` are wrapper functions of
[`fastcpd()`](https://fastcpd.xingchi.li/reference/fastcpd.md) to find
change points in ARIMA(\\p\\, \\d\\, \\q\\) models. The function is
similar to
[`fastcpd()`](https://fastcpd.xingchi.li/reference/fastcpd.md) except
that the data is by default a one-column matrix or univariate vector and
thus a formula is not required here.

## Usage

``` r
fastcpd_arima(data, order = 0, ...)

fastcpd.arima(data, order = 0, ...)
```

## Arguments

- data:

  A numeric vector, a matrix, a data frame or a time series object.

- order:

  A vector of length three specifying the order of the ARIMA model.

- ...:

  Other arguments passed to
  [`fastcpd()`](https://fastcpd.xingchi.li/reference/fastcpd.md), for
  example, `segment_count`. One special argument can be passed here is
  `include.mean`, which is a logical value indicating whether the mean
  should be included in the model. The default value is `TRUE`.

## Value

A [fastcpd](https://fastcpd.xingchi.li/reference/fastcpd-class.md)
object.

## See also

[`fastcpd()`](https://fastcpd.xingchi.li/reference/fastcpd.md)

## Examples

``` r
# \donttest{
set.seed(1)
n <- 801
w <- rnorm(n + 1, 0, 3)
dx <- rep(0, n + 1)
x <- rep(0, n + 1)
for (i in 1:400) {
  dx[i + 1] <- 0.9 * dx[i] + w[i + 1] - 0.1 * w[i]
  x[i + 1] <- x[i] + dx[i + 1]
}
for (i in 401:n) {
  dx[i + 1] <- -0.6 * dx[i] + w[i + 1] + 0.3 * w[i]
  x[i + 1] <- x[i] + dx[i + 1]
}
result <- fastcpd.arima(
  diff(x[1 + seq_len(n)]),
  c(1, 0, 1),
  segment_count = 3,
  include.mean = FALSE
)
#> Warning: possible convergence problem: optim gave code = 1
#> Warning: possible convergence problem: optim gave code = 1
#> Warning: possible convergence problem: optim gave code = 1
#> Warning: possible convergence problem: optim gave code = 1
#> Warning: possible convergence problem: optim gave code = 1
#> Warning: possible convergence problem: optim gave code = 1
#> Warning: possible convergence problem: optim gave code = 1
#> Warning: possible convergence problem: optim gave code = 1
#> Warning: possible convergence problem: optim gave code = 1
#> Warning: possible convergence problem: optim gave code = 1
#> Warning: possible convergence problem: optim gave code = 1
#> Warning: possible convergence problem: optim gave code = 1
#> Warning: possible convergence problem: optim gave code = 1
#> Warning: possible convergence problem: optim gave code = 1
#> Warning: possible convergence problem: optim gave code = 1
#> Warning: possible convergence problem: optim gave code = 1
#> Warning: possible convergence problem: optim gave code = 1
summary(result)
#> 
#> Call:
#> fastcpd.arima(data = diff(x[1 + seq_len(n)]), order = c(1, 0, 
#>     1), segment_count = 3, include.mean = FALSE)
#> 
#> Change points:
#> 429 
#> 
#> Cost values:
#> 1088.875 973.5151 
#> 
#> Parameters:
#>   segment 1 segment 2
#> 1         0         0
#> 2         0         0
#> 3         0         0
plot(result)

# }
```
