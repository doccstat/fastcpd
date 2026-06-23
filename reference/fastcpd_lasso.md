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
#> 191.8539 245.4502 242.25 282.737 
#> 
#> Parameters:
#>       segment 1   segment 2   segment 3     segment 4
#> 1  -2.473168539 -0.98447207  3.46244792 -2.4626421537
#> 2  -2.588593610 -2.41828115  2.10848311  3.7541930968
#> 3  -3.983941649 -2.18197691  4.59036172 -1.0198868501
#> 4  -4.378913735  0.91071392  3.95765503 -2.8685363860
#> 5  -4.495473996  0.00000000  4.84060707  2.7550202120
#> 6   0.000000000  0.00000000  0.00000000  0.0000000000
#> 7  -0.053395310  0.00000000  0.00000000  0.0000000000
#> 8   0.000000000  0.00000000  0.00000000  0.0000000000
#> 9   0.000000000  0.00000000  0.00000000  0.0000000000
#> 10  0.000000000  0.00000000  0.00000000  0.0000000000
#> 11  0.000000000  0.00000000  0.00000000  0.0000000000
#> 12  0.000000000  0.00000000  0.00000000  0.0000000000
#> 13  0.000000000  0.00000000  0.00000000  0.0000000000
#> 14  0.000000000  0.00000000  0.00000000  0.0000000000
#> 15  0.000000000  0.00000000  0.00000000  0.0000000000
#> 16  0.000000000  0.00000000  0.00000000  0.0000000000
#> 17  0.000000000  0.00000000  0.00000000  0.0000000000
#> 18  0.000000000  0.00000000  0.00000000  0.0000000000
#> 19  0.000000000  0.00000000  0.00000000  0.0000000000
#> 20  0.000000000  0.00000000  0.00000000  0.0000000000
#> 21  0.000000000  0.00000000  0.00000000  0.0000000000
#> 22  0.015707059  0.00000000 -0.09845551  0.0000000000
#> 23  0.000000000  0.00000000  0.00000000  0.0000000000
#> 24  0.000000000  0.00000000  0.00000000  0.0000000000
#> 25  0.000000000  0.00000000  0.00000000  0.0000000000
#> 26  0.000000000 -0.09900977  0.00000000  0.0000000000
#> 27  0.000000000  0.00000000  0.00000000  0.0000000000
#> 28  0.000000000  0.00000000  0.00000000  0.0000000000
#> 29  0.000000000  0.00000000  0.00000000  0.0000000000
#> 30  0.000000000  0.00000000  0.00000000  0.0000000000
#> 31  0.000000000  0.00000000  0.00000000  0.0000000000
#> 32  0.000000000  0.00000000  0.00000000  0.0000000000
#> 33  0.177142836  0.00000000  0.00000000  0.0000000000
#> 34  0.000000000  0.00000000  0.00000000  0.0000000000
#> 35  0.000000000  0.00000000  0.00000000  0.0000000000
#> 36  0.005010027  0.00000000  0.00000000 -0.0002328763
#> 37  0.000000000  0.00000000  0.00000000  0.0000000000
#> 38  0.000000000  0.00000000  0.00000000  0.0000000000
#> 39  0.000000000  0.00000000  0.00000000  0.0000000000
#> 40  0.000000000 -0.08092709  0.00000000  0.0000000000
#> 41  0.000000000  0.00000000  0.00000000  0.0000000000
#> 42  0.000000000  0.00000000  0.00000000  0.0000000000
#> 43  0.000000000  0.00000000  0.00000000  0.0000000000
#> 44  0.000000000  0.00000000  0.00000000  0.0000000000
#> 45 -0.037051417  0.00000000  0.00000000  0.0000000000
#> 46  0.000000000  0.00000000  0.00000000  0.0000000000
#> 47  0.000000000  0.00000000  0.00000000  0.0000000000
#> 48  0.000000000  0.00000000  0.00000000  0.0000000000
#> 49  0.000000000  0.00000000  0.00000000  0.0000000000
#> 50  0.000000000  0.00000000  0.00000000  0.0000000000


# }
```
