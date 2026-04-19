# Find change points efficiently in Poisson regression models

`fastcpd_poisson()` and `fastcpd.poisson()` are wrapper functions of
[`fastcpd()`](https://fastcpd.xingchi.li/reference/fastcpd.md) to find
change points in Poisson regression models. The function is similar to
[`fastcpd()`](https://fastcpd.xingchi.li/reference/fastcpd.md) except
that the data is by default a matrix or data frame with the response
variable as the first column and thus a formula is not required here.

## Usage

``` r
fastcpd_poisson(data, ...)

fastcpd.poisson(data, ...)
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
# \donttest{
if (requireNamespace("mvtnorm", quietly = TRUE)) {
  set.seed(1)
  n <- 1100
  p <- 3
  x <- mvtnorm::rmvnorm(n, rep(0, p), diag(p))
  delta <- rnorm(p)
  theta_0 <- c(1, 0.3, -1)
  y <- c(
    rpois(500, exp(x[1:500, ] %*% theta_0)),
    rpois(300, exp(x[501:800, ] %*% (theta_0 + delta))),
    rpois(200, exp(x[801:1000, ] %*% theta_0)),
    rpois(100, exp(x[1001:1100, ] %*% (theta_0 - delta)))
  )
  result <- fastcpd.poisson(cbind(y, x))
  summary(result)
  plot(result)
}
#> 
#> Call:
#> fastcpd.poisson(data = cbind(y, x))
#> 
#> Change points:
#> 506 838 1003 
#> 
#> Cost values:
#> 248.1639 221.1317 74.88469 48.95659 
#> 
#> Parameters:
#>    segment 1  segment 2  segment 3  segment 4
#> 1  1.0154681  0.6568705  1.0371861  1.4451928
#> 2  0.2763783 -0.2131387  0.2648813  0.9910079
#> 3 -1.0493262 -0.5942795 -0.9801554 -1.4354638

# }
```
