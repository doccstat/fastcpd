# Find change points efficiently in logistic regression models

`fastcpd_binomial()` and `fastcpd.binomial()` are wrapper functions of
[`fastcpd()`](https://fastcpd.xingchi.li/reference/fastcpd.md) to find
change points in logistic regression models. The function is similar to
[`fastcpd()`](https://fastcpd.xingchi.li/reference/fastcpd.md) except
that the data is by default a matrix or data frame with the response
variable as the first column and thus a formula is not required here.

## Usage

``` r
fastcpd_binomial(data, ...)

fastcpd.binomial(data, ...)
```

## Arguments

- data:

  A matrix or a data frame with the response variable as the first
  column.

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
if (requireNamespace("mvtnorm", quietly = TRUE)) {
  set.seed(1)
  n <- 500
  p <- 4
  x <- mvtnorm::rmvnorm(n, rep(0, p), diag(p))
  theta <- rbind(rnorm(p, 0, 1), rnorm(p, 2, 1))
  y <- c(
    rbinom(300, 1, 1 / (1 + exp(-x[1:300, ] %*% theta[1, ]))),
    rbinom(200, 1, 1 / (1 + exp(-x[301:n, ] %*% theta[2, ])))
  )
  result <- suppressWarnings(fastcpd.binomial(cbind(y, x)))
  summary(result)
  plot(result)
}
#> 
#> Call:
#> fastcpd.binomial(data = cbind(y, x))
#> 
#> Change points:
#> 302 
#> 
#> Cost values:
#> 136.8846 66.69302 
#> 
#> Parameters:
#>    segment 1 segment 2
#> 1 -0.9260182 2.1294962
#> 2 -1.6033835 2.7583247
#> 3  1.0343338 2.3818010
#> 4  0.3653870 0.7261152
```
