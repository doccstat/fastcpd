# Find change points efficiently in mean variance change models

`fastcpd_meanvariance()`, `fastcpd.meanvariance()`, `fastcpd_mv()`,
`fastcpd.mv()` are wrapper functions of
[`fastcpd()`](https://fastcpd.xingchi.li/reference/fastcpd.md) to find
the meanvariance change. The function is similar to
[`fastcpd()`](https://fastcpd.xingchi.li/reference/fastcpd.md) except
that the data is by default a matrix or data frame or a vector with each
row / element as an observation and thus a formula is not required here.

## Usage

``` r
fastcpd_meanvariance(data, ...)

fastcpd.meanvariance(data, ...)

fastcpd_mv(data, ...)

fastcpd.mv(data, ...)
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
data <- c(
  rnorm(3000, 0, 1),
  rnorm(1000, 10, 1),
  rnorm(3000, 10, 20),
  rnorm(1000, 0, 1)
)
system.time(result <- fastcpd.mv(data))
#>    user  system elapsed 
#>   0.006   0.002   0.008 
summary(result)
#> 
#> Call:
#> fastcpd.mv(data = data)
#> 
#> Change points:
#> 3000 4000 7000 
#> 
#> Cost values:
#> 110.9701 44.33282 9001.941 7.441331 
#> 
#> Parameters:
#>     segment 1 segment 2  segment 3    segment 4
#> 1 -0.00420034 10.016722   9.650309 -0.004211938
#> 2  1.07141035  1.078801 401.934668  1.002069791
plot(result)

set.seed(1)
p <- 3
data <- if (requireNamespace("mvtnorm", quietly = TRUE)) {
  rbind(
    mvtnorm::rmvnorm(2e+5, mean = rep(0, p), sigma = diag(1, p)),
    mvtnorm::rmvnorm(1e+5, mean = rep(50, p), sigma = diag(1, p)),
    mvtnorm::rmvnorm(2e+5, mean = rep(0, p), sigma = diag(100, p)),
    mvtnorm::rmvnorm(2e+5, mean = rep(0, p), sigma = diag(1, p)),
    mvtnorm::rmvnorm(1e+5, mean = rep(50, p), sigma = diag(1, p)),
    mvtnorm::rmvnorm(2e+5, mean = rep(50, p), sigma = diag(100, p))
  )
} else {
  rbind(
    matrix(rnorm(p * 2e+5, mean = 0, sd = 1), ncol = p),
    matrix(rnorm(p * 1e+5, mean = 50, sd = 1), ncol = p),
    matrix(rnorm(p * 2e+5, mean = 0, sd = 10), ncol = p),
    matrix(rnorm(p * 2e+5, mean = 0, sd = 1), ncol = p),
    matrix(rnorm(p * 1e+5, mean = 50, sd = 1), ncol = p),
    matrix(rnorm(p * 2e+5, mean = 50, sd = 10), ncol = p)
  )
}
system.time(result <- fastcpd.mv(data))
#>    user  system elapsed 
#>   1.314   0.064   1.378 
summary(result)
#> 
#> Call:
#> fastcpd.mv(data = data)
#> 
#> Change points:
#> 2e+05 3e+05 499999 7e+05 8e+05 
#> 
#> Cost values:
#> 108.8674 -43.02536 1382904 -47.68093 377.4763 1381660 
#> 
#> Parameters:
#>       segment 1     segment 2     segment 3     segment 4    segment 5
#> 1   0.001011082 50.0003528125 -1.647663e-02  0.0032725571 49.995348620
#> 2  -0.001572291 50.0047380197  2.735846e-02 -0.0002665615 49.995608790
#> 3  -0.001048280 49.9949663194 -2.119325e-04 -0.0023005354 50.004925725
#> 4   1.001716184  1.0053160771  1.000400e+02  1.0049799581  1.011184355
#> 5   0.003146441 -0.0061509213  1.558501e-01  0.0031928989  0.002837378
#> 6   0.001702982 -0.0009432632 -1.746045e-01 -0.0030281680 -0.001092170
#> 7   0.003146441 -0.0061509213  1.558501e-01  0.0031928989  0.002837378
#> 8   0.998570871  0.9954591727  1.004605e+02  0.9942135351  0.999139034
#> 9   0.001008951 -0.0062795529  7.083125e-02  0.0006677883 -0.001649613
#> 10  0.001702982 -0.0009432632 -1.746045e-01 -0.0030281680 -0.001092170
#> 11  0.001008951 -0.0062795529  7.083125e-02  0.0006677883 -0.001649613
#> 12  1.000100554  0.9971197656  1.007926e+02  0.9996613540  0.995957111
#>       segment 6
#> 1   50.00765767
#> 2   49.97689987
#> 3   50.01812646
#> 4  100.43738403
#> 5    0.09929858
#> 6    0.04291338
#> 7    0.09929858
#> 8   99.71858356
#> 9   -0.09614237
#> 10   0.04291338
#> 11  -0.09614237
#> 12  99.88324093
result@thetas[seq_len(p), ]
#>      segment 1 segment 2     segment 3     segment 4 segment 5 segment 6
#> 1  0.001011082  50.00035 -0.0164766285  0.0032725571  49.99535  50.00766
#> 2 -0.001572291  50.00474  0.0273584578 -0.0002665615  49.99561  49.97690
#> 3 -0.001048280  49.99497 -0.0002119325 -0.0023005354  50.00493  50.01813
lapply(result@thetas[seq_len(p^2) + p, ], function(thetas) matrix(thetas, p))
#> $`segment 1`
#>             [,1]        [,2]        [,3]
#> [1,] 1.001716184 0.003146441 0.001702982
#> [2,] 0.003146441 0.998570871 0.001008951
#> [3,] 0.001702982 0.001008951 1.000100554
#> 
#> $`segment 2`
#>               [,1]         [,2]          [,3]
#> [1,]  1.0053160771 -0.006150921 -0.0009432632
#> [2,] -0.0061509213  0.995459173 -0.0062795529
#> [3,] -0.0009432632 -0.006279553  0.9971197656
#> 
#> $`segment 3`
#>             [,1]         [,2]         [,3]
#> [1,] 100.0399627   0.15585014  -0.17460453
#> [2,]   0.1558501 100.46054935   0.07083125
#> [3,]  -0.1746045   0.07083125 100.79256474
#> 
#> $`segment 4`
#>              [,1]         [,2]          [,3]
#> [1,]  1.004979958 0.0031928989 -0.0030281680
#> [2,]  0.003192899 0.9942135351  0.0006677883
#> [3,] -0.003028168 0.0006677883  0.9996613540
#> 
#> $`segment 5`
#>              [,1]         [,2]         [,3]
#> [1,]  1.011184355  0.002837378 -0.001092170
#> [2,]  0.002837378  0.999139034 -0.001649613
#> [3,] -0.001092170 -0.001649613  0.995957111
#> 
#> $`segment 6`
#>              [,1]        [,2]        [,3]
#> [1,] 100.43738403  0.09929858  0.04291338
#> [2,]   0.09929858 99.71858356 -0.09614237
#> [3,]   0.04291338 -0.09614237 99.88324093
#> 
```
