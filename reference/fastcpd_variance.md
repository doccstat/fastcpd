# Find change points efficiently in variance change models

`fastcpd_variance()` and `fastcpd.variance()` are wrapper functions of
[`fastcpd()`](https://fastcpd.xingchi.li/reference/fastcpd.md) to find
the variance change. The function is similar to
[`fastcpd()`](https://fastcpd.xingchi.li/reference/fastcpd.md) except
that the data is by default a matrix or data frame or a vector with each
row / element as an observation and thus a formula is not required here.

## Usage

``` r
fastcpd_variance(data, ...)

fastcpd.variance(data, ...)
```

## Arguments

- data:

  A matrix, a data frame or a vector.

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
data <- c(rnorm(300, 0, 1), rnorm(400, 0, 10), rnorm(300, 0, 1))
system.time(result <- fastcpd.variance(data))
#>    user  system elapsed 
#>   0.003   0.000   0.004 
summary(result)
#> 
#> Call:
#> fastcpd.variance(data = data)
#> 
#> Change points:
#> 300 700 
#> 
#> Cost values:
#> 1.341031 945.1872 30.94909 
#> 
#> Parameters:
#>   segment 1 segment 2 segment 3
#> 1 0.9287098  111.3079  1.159284
plot(result)

if (requireNamespace("mvtnorm", quietly = TRUE)) {
  set.seed(1)
  p <- 3
  data <- rbind(
    mvtnorm::rmvnorm(
      3e+5, rep(0, p), crossprod(matrix(runif(p^2) * 2 - 1, p))
    ),
    mvtnorm::rmvnorm(
      4e+5, rep(0, p), crossprod(matrix(runif(p^2) * 2 - 1, p))
    ),
    mvtnorm::rmvnorm(
      3e+5, rep(0, p), crossprod(matrix(runif(p^2) * 2 - 1, p))
    )
  )
  result_time <- system.time(
    result <- fastcpd.variance(data, r.progress = FALSE, cp_only = TRUE)
  )
  print(result_time)
  summary(result)
}
#>    user  system elapsed 
#>   0.671   0.053   0.724 
#> 
#> Call:
#> fastcpd.variance(data = data, r.progress = FALSE, cp_only = TRUE)
#> 
#> Change points:
#> 300002 700003 
```
