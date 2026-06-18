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
#> 178.2183 238.3949 242.25 297.784 
#> 
#> Parameters:
#>      segment 1  segment 2   segment 3 segment 4
#> 1  -2.51973417 -1.0190850  3.46244792 -2.405264
#> 2  -2.63509034 -2.4531499  2.10848311  3.712577
#> 3  -4.05586278 -2.2110852  4.59036172 -0.961646
#> 4  -4.37980262  0.9340838  3.95765503 -2.809584
#> 5  -4.55278268  0.0000000  4.84060707  2.715393
#> 6   0.00000000  0.0000000  0.00000000  0.000000
#> 7  -0.07038421  0.0000000  0.00000000  0.000000
#> 8   0.00000000  0.0000000  0.00000000  0.000000
#> 9   0.00000000  0.0000000  0.00000000  0.000000
#> 10  0.00000000  0.0000000  0.00000000  0.000000
#> 11  0.00000000  0.0000000  0.00000000  0.000000
#> 12  0.00000000  0.0000000  0.00000000  0.000000
#> 13  0.00000000  0.0000000  0.00000000  0.000000
#> 14  0.00000000  0.0000000  0.00000000  0.000000
#> 15  0.00000000  0.0000000  0.00000000  0.000000
#> 16  0.00000000  0.0000000  0.00000000  0.000000
#> 17  0.00000000  0.0000000  0.00000000  0.000000
#> 18  0.00000000  0.0000000  0.00000000  0.000000
#> 19  0.00000000  0.0000000  0.00000000  0.000000
#> 20  0.02787591  0.0000000  0.00000000  0.000000
#> 21  0.00000000  0.0000000  0.00000000  0.000000
#> 22  0.06169046  0.0000000 -0.09845551  0.000000
#> 23  0.00000000  0.0000000  0.00000000  0.000000
#> 24  0.00000000  0.0000000  0.00000000  0.000000
#> 25  0.00000000  0.0000000  0.00000000  0.000000
#> 26  0.00000000 -0.1146819  0.00000000  0.000000
#> 27  0.00000000  0.0000000  0.00000000  0.000000
#> 28  0.00000000  0.0000000  0.00000000  0.000000
#> 29  0.00000000  0.0000000  0.00000000  0.000000
#> 30  0.00000000  0.0000000  0.00000000  0.000000
#> 31  0.00000000  0.0000000  0.00000000  0.000000
#> 32  0.00000000  0.0000000  0.00000000  0.000000
#> 33  0.20636226  0.0000000  0.00000000  0.000000
#> 34  0.01366932  0.0000000  0.00000000  0.000000
#> 35  0.00000000  0.0000000  0.00000000  0.000000
#> 36  0.06957124  0.0000000  0.00000000  0.000000
#> 37  0.00000000  0.0000000  0.00000000  0.000000
#> 38  0.00000000  0.0000000  0.00000000  0.000000
#> 39  0.00000000  0.0000000  0.00000000  0.000000
#> 40  0.00000000 -0.1038122  0.00000000  0.000000
#> 41  0.00000000  0.0000000  0.00000000  0.000000
#> 42  0.00000000  0.0000000  0.00000000  0.000000
#> 43  0.00000000  0.0000000  0.00000000  0.000000
#> 44  0.00000000  0.0000000  0.00000000  0.000000
#> 45 -0.06216254  0.0000000  0.00000000  0.000000
#> 46  0.00000000  0.0000000  0.00000000  0.000000
#> 47  0.00000000  0.0000000  0.00000000  0.000000
#> 48 -0.03308599  0.0000000  0.00000000  0.000000
#> 49  0.03360350  0.0000000  0.00000000  0.000000
#> 50  0.00000000  0.0000000  0.00000000  0.000000


# }
```
