# Find change points efficiently in VAR(\\p\\) models

`fastcpd_var()` and `fastcpd.var()` are wrapper functions of
[`fastcpd_ts()`](https://fastcpd.xingchi.li/reference/fastcpd_ts.md) to
find change points in VAR(\\p\\) models. The function is similar to
[`fastcpd_ts()`](https://fastcpd.xingchi.li/reference/fastcpd_ts.md)
except that the data is by default a matrix with row as an observation
and thus a formula is not required here.

## Usage

``` r
fastcpd_var(data, order = 0, ...)

fastcpd.var(data, order = 0, ...)
```

## Arguments

- data:

  A matrix, a data frame or a time series object.

- order:

  A positive integer specifying the order of the VAR model.

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
n <- 300
p <- 2
theta_1 <- matrix(c(-0.3, 0.6, -0.5, 0.4, 0.2, 0.2, 0.2, -0.2), nrow = p)
theta_2 <- matrix(c(0.3, -0.4, 0.1, -0.5, -0.5, -0.2, -0.5, 0.2), nrow = p)
x <- matrix(0, n + 2, p)
for (i in 1:200) {
  x[i + 2, ] <- theta_1 %*% c(x[i + 1, ], x[i, ]) + rnorm(p, 0, 1)
}
for (i in 201:n) {
  x[i + 2, ] <- theta_2 %*% c(x[i + 1, ], x[i, ]) + rnorm(p, 0, 1)
}
result <- fastcpd.var(x, 2)
summary(result)
#> 
#> Call:
#> fastcpd.var(data = x, order = 2)
#> 
#> Change points:
#> 204 
#> 
#> Cost values:
#> 583.9159 306.6093 
#> 
#> Parameters:
#>    segment 1  segment 2
#> 1 -0.2524511  0.2760795
#> 2 -0.5679217  0.2764738
#> 3  0.2065672 -0.5327953
#> 4  0.2451136 -0.4303701
#> 5  0.5655444 -0.3844917
#> 6  0.3903103 -0.5302748
#> 7  0.1614109 -0.1918981
#> 8 -0.1716114  0.1671678
```
