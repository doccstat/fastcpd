# Find change points efficiently in penalized linear regression models

`fastcpd_lasso()` and `fastcpd.lasso()` are wrapper functions of
[`fastcpd()`](https://fastcpd.xingchi.li/reference/fastcpd.md) to find
change points in penalized linear regression models. The function is
similar to
[`fastcpd()`](https://fastcpd.xingchi.li/reference/fastcpd.md) except
that the data is by default a matrix or data frame with the response
variable as the first column and thus a formula is not required here.

## Usage

``` r
fastcpd_lasso(data, ...)

fastcpd.lasso(data, ...)
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
if (
  requireNamespace("ggplot2", quietly = TRUE) &&
    requireNamespace("mvtnorm", quietly = TRUE)
) {
  set.seed(1)
  n <- 480
  p_true <- 5
  p <- 50
  x <- mvtnorm::rmvnorm(n, rep(0, p), diag(p))
  theta_0 <- rbind(
    runif(p_true, -5, -2),
    runif(p_true, -3, 3),
    runif(p_true, 2, 5),
    runif(p_true, -5, 5)
  )
  theta_0 <- cbind(theta_0, matrix(0, ncol = p - p_true, nrow = 4))
  y <- c(
    x[1:80, ] %*% theta_0[1, ] + rnorm(80, 0, 1),
    x[81:200, ] %*% theta_0[2, ] + rnorm(120, 0, 1),
    x[201:320, ] %*% theta_0[3, ] + rnorm(120, 0, 1),
    x[321:n, ] %*% theta_0[4, ] + rnorm(160, 0, 1)
  )
  result <- fastcpd.lasso(
    cbind(y, x),
    multiple_epochs = function(segment_length) if (segment_length < 30) 1 else 0
  )
  summary(result)
  plot(result)

  # Combine estimated thetas with true parameters
  thetas <- result@thetas
  thetas <- cbind.data.frame(thetas, t(theta_0))
  names(thetas) <- c(
    "segment 1", "segment 2", "segment 3", "segment 4",
    "segment 1 truth", "segment 2 truth", "segment 3 truth", "segment 4 truth"
  )
  thetas$coordinate <- c(seq_len(p_true), rep("rest", p - p_true))

  # Melt the data frame using base R (i.e., convert from wide to long format)
  data_cols <- setdiff(names(thetas), "coordinate")
  molten <- data.frame(
    coordinate = rep(thetas$coordinate, times = length(data_cols)),
    variable = rep(data_cols, each = nrow(thetas)),
    value = as.vector(as.matrix(thetas[, data_cols]))
  )

  # Remove the "segment " and " truth" parts to extract the segment number
  molten$segment <- gsub("segment ", "", molten$variable)
  molten$segment <- gsub(" truth", "", molten$segment)

  # Compute height: the numeric value of the segment plus an offset for truth values
  molten$height <- as.numeric(gsub("segment.*", "", molten$segment)) +
    0.2 * as.numeric(grepl("truth", molten$variable))

  # Create a parameter indicator based on whether the variable corresponds to truth or estimation
  molten$parameter <- ifelse(grepl("truth", molten$variable), "truth", "estimated")

  p <- ggplot2::ggplot() +
    ggplot2::geom_point(
      data = molten,
      ggplot2::aes(
        x = value, y = height, shape = coordinate, color = parameter
      ),
      size = 4
    ) +
    ggplot2::ylim(0.8, 4.4) +
    ggplot2::ylab("segment") +
    ggplot2::theme_bw()
  print(p)
}
#> 
#> Call:
#> fastcpd.lasso(data = cbind(y, x), multiple_epochs = function(segment_length) if (segment_length < 
#>     30) 1 else 0)
#> 
#> Change points:
#> 80 200 321 
#> 
#> Cost values:
#> 204.8206 244.7625 244.5833 284.7771 
#> 
#> Parameters:
#> 50 x 4 sparse Matrix of class "dgCMatrix"
#>          segment 1   segment 2   segment 3    segment 4
#>  [1,] -2.430099838 -0.98559186  3.44703078 -2.447825206
#>  [2,] -2.558517191 -2.42611597  2.09094273  3.736499163
#>  [3,] -3.914451531 -2.17771957  4.59365835 -1.014029107
#>  [4,] -4.327584111  0.90347628  3.94608206 -2.871825706
#>  [5,] -4.444653183  .           4.84648763  2.738991414
#>  [6,]  .            .           .           0.003729822
#>  [7,]  .            .           .           .          
#>  [8,]  .            .           .           .          
#>  [9,]  .            .           .           .          
#> [10,]  .            .           .           .          
#> [11,]  .            .           .           .          
#> [12,]  .            .           .           .          
#> [13,]  .            .           .           .          
#> [14,]  .            .           .           .          
#> [15,]  .            .           .           .          
#> [16,]  .            .           .           .          
#> [17,]  .            .           .           .          
#> [18,]  .            .           .           .          
#> [19,]  .            .           .           .          
#> [20,]  .            .           .           .          
#> [21,]  .            .           .           .          
#> [22,]  0.005165931  .          -0.07024678  .          
#> [23,]  .            .           .           .          
#> [24,]  .            .           .           .          
#> [25,]  .            .           .           .          
#> [26,]  .           -0.10792142  .           .          
#> [27,]  .            .           .           .          
#> [28,]  .            .           .           .          
#> [29,]  .            .           .           .          
#> [30,]  .            .           .           .          
#> [31,]  .            .           .           .          
#> [32,]  .            .           .           .          
#> [33,]  0.114429784  .           .           .          
#> [34,]  .            .           .           .          
#> [35,]  .            .           .           .          
#> [36,]  .            .           .          -0.002930384
#> [37,]  .            .           .           .          
#> [38,]  .            .           .           .          
#> [39,]  .            .           .           .          
#> [40,]  .           -0.08532151  .           .          
#> [41,]  .            .           .           .          
#> [42,]  .            .           .           .          
#> [43,]  .            .           .           .          
#> [44,]  .            .           .           .          
#> [45,]  .            .           .           .          
#> [46,]  .            .           .           .          
#> [47,]  .            .           .           .          
#> [48,]  .            .           .           .          
#> [49,]  .            .           .           .          
#> [50,]  .            .           .           .          


# }
```
