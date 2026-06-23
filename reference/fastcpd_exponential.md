# Find change points efficiently in exponentially distributed data

`fastcpd_exponential()` and `fastcpd.exponential()` are wrapper
functions of
[`fastcpd()`](https://fastcpd.xingchi.li/reference/fastcpd.md) to find
changes in the rate of exponentially distributed data, i.e. mean change
under exponentially distributed noise (cf. `changepoint::cpt.meanvar`
with `test.stat = "Exponential"`). The function is similar to
[`fastcpd()`](https://fastcpd.xingchi.li/reference/fastcpd.md) except
that the data is by default a matrix or data frame or a vector with each
row / element as an observation and thus a formula is not required here.

## Usage

``` r
fastcpd_exponential(data, ...)

fastcpd.exponential(data, ...)
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
  rexp(300, rate = 1),
  rexp(400, rate = 0.1),
  rexp(300, rate = 1)
))
system.time(result <- fastcpd.exponential(data))
#>    user  system elapsed 
#>   0.003   0.000   0.002 
summary(result)
#> 
#> Call:
#> fastcpd.exponential(data = data)
#> 
#> Change points:
#> 300 700 
#> 
#> Cost values:
#> 301.6378 1326.448 331.5256 
#> 
#> Parameters:
#>   segment 1  segment 2 segment 3
#> 1  1.004055 0.09939726 0.9088465
plot(result)
```
