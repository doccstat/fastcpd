# Find change points efficiently via rank-based change point detection

`fastcpd_rank()` and `fastcpd.rank()` detect change points using a
rank-based, distribution-free cost. Each column is replaced by its
global rank centred at zero, and change points in the mean of these
centred ranks are detected with the existing PELT infrastructure. The
result is fully deterministic and requires no bandwidth selection. The
method is most powerful for location shifts; for scale-only or general
distributional changes,
[`fastcpd_kcp()`](https://fastcpd.xingchi.li/reference/fastcpd_kcp.md)
is preferable.

## Usage

``` r
fastcpd_rank(data, ...)

fastcpd.rank(data, ...)
```

## Arguments

- data:

  A numeric vector or a matrix with one row per observation.

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
x <- c(rt(200, df = 2), rt(200, df = 2, ncp = 5))
result_mean <- fastcpd_mean(x)
summary(result_mean)
#> 
#> Call:
#> fastcpd_mean(data = x)
#> 
#> Change points:
#> 332 333 
#> 
#> Cost values:
#> 3.352105 0 2.150836 
#> 
#> Parameters:
#>   segment 1 segment 2 segment 3
#> 1  3.436182  3047.762  8.251482
result_rank <- fastcpd_rank(x)
summary(result_rank)
#> 
#> Call:
#> fastcpd_rank(data = x)
#> 
#> Change points:
#> 200 
#> 
#> Cost values:
#> 102.8266 85.86544 
#> 
#> Parameters:
#>   segment 1 segment 2
#> 1    -94.44     94.44
```
