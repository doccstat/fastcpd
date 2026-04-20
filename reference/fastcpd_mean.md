# Find change points efficiently in mean change models

`fastcpd_mean()` and `fastcpd.mean()` are wrapper functions of
[`fastcpd()`](https://fastcpd.xingchi.li/reference/fastcpd.md) to find
the mean change. The function is similar to
[`fastcpd()`](https://fastcpd.xingchi.li/reference/fastcpd.md) except
that the data is by default a matrix or data frame or a vector with each
row / element as an observation and thus a formula is not required here.

## Usage

``` r
fastcpd_mean(data, ...)

fastcpd.mean(data, ...)
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
data <- matrix(c(
  rnorm(300, mean = 0, sd = 10),
  rnorm(400, mean = 50, sd = 10),
  rnorm(300, mean = 2, sd = 10)
))
system.time(result <- fastcpd.mean(data))
#>    user  system elapsed 
#>   0.006   0.000   0.005 
summary(result)
#> 
#> Call:
#> fastcpd.mean(data = data)
#> 
#> Change points:
#> 300 700 
#> 
#> Cost values:
#> 124.2113 197.0938 154.3416 
#> 
#> Parameters:
#>   segment 1 segment 2 segment 3
#> 1 0.3358428  49.42092  2.047991
plot(result)

set.seed(1)
p <- 3
data <- rbind(
  matrix(rnorm(p * 3e+5, mean = 0, sd = 10), ncol = p),
  matrix(rnorm(p * 4e+5, mean = 50, sd = 10), ncol = p),
  matrix(rnorm(p * 3e+5, mean = 2, sd = 10), ncol = p)
)
system.time(result <- fastcpd.mean(data, r.progress = FALSE, cp_only = TRUE))
#>    user  system elapsed 
#>   3.852   0.310   3.856 
summary(result)
#> 
#> Call:
#> fastcpd.mean(data = data, r.progress = FALSE, cp_only = TRUE)
#> 
#> Change points:
#> 3e+05 7e+05 
```
