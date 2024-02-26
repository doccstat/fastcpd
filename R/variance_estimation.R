#' @title Variance estimation for linear models with change points
#'
#' @description Estimate the variance for each block and then take the average.
#'
#' @example tests/testthat/examples/variance_lm.R
#'
#' @md
#'
#' @param data A matrix or a data frame with the response variable as the first
#'   column.
#' @param block_size The size of the blocks to use for variance estimation.
#' @param p The number of predictors in the model.
#'
#' @return A numeric value representing the variance.
#'
#' @rdname variance_lm
#' @export
variance_lm <- function(data, block_size = NULL, p = NULL) {
  if (is.null(p)) {
    p <- ncol(data) - 1
  }
  if (is.null(block_size)) {
    block_size <- p + 1
  }
  n <- nrow(data)
  variance_estimation <- rep(NA, n - block_size)
  for (i in seq_len(n - block_size)) {
    block_index <- seq_len(block_size) + i - 1
    block_index_lagged <- seq_len(block_size) + i

    y_block <- data[block_index, 1]
    x_block <- data[block_index, -1, drop = FALSE]

    y_block_lagged <- data[block_index_lagged, 1]
    x_block_lagged <- data[block_index_lagged, -1, drop = FALSE]

    x_t_x <- crossprod(x_block)
    x_t_x_lagged <- crossprod(x_block_lagged)

    tryCatch(
      expr = {
        block_slope <-
          solve(crossprod(x_block), crossprod(x_block, y_block))
        block_lagged_slope <- solve(
          crossprod(x_block_lagged),
          crossprod(x_block_lagged, y_block_lagged)
        )
        x_t_x_inv <- solve(x_t_x)
        x_t_x_inv_lagged <- solve(x_t_x_lagged)
        inv_product <- x_t_x_inv %*% x_t_x_inv_lagged
        cross_term_x <- crossprod(
          x_block[-1, , drop = FALSE],
          x_block_lagged[-block_size, , drop = FALSE]
        )
        cross_term <- inv_product %*% cross_term_x
        delta_numerator <- crossprod(block_slope - block_lagged_slope)
        delta_denominator <-
          sum(diag(x_t_x_inv + x_t_x_inv_lagged - 2 * cross_term))
        variance_estimation[i] <- delta_numerator / delta_denominator
      },
      error = function(e) {
        variance_estimation[i] <- NA
        message("Variance estimation failed for block ", i, ".")
      }
    )
  }
  mean(variance_estimation, na.rm = TRUE)
}

#' @rdname variance_lm
#' @export
variance.lm <- variance_lm  # nolint: Conventional R function style

#' @title Variance estimation for mean change models
#'
#' @description Implement Rice estimator.
#'
#' @example tests/testthat/examples/variance_mean.R
#'
#' @md
#'
#' @param data A matrix or a data frame with data points as each row.
#'
#' @return A matrix representing the variance-covariance matrix or a numeric
#'   value representing the variance.
#'
#' @rdname variance_mean
#' @export
variance_mean <- function(data) {
  n <- nrow(data)
  p <- ncol(data)
  variance_estimation <- array(NA, c(n - 1, p, p))
  for (i in seq_len(n - 1)) {
    variance_estimation[i, , ] <- crossprod(
      data[i + 1, , drop = FALSE] - data[i, , drop = FALSE]
    )
  }
  sigma <- colMeans(variance_estimation, na.rm = TRUE) / 2
  if (p == 1) {
    c(sigma)
  } else {
    sigma
  }
}

#' @rdname variance_mean
#' @export
variance.mean <- variance_mean  # nolint: Conventional R function style
