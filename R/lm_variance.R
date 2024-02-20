#' @title Variance estimation for linear models with change points
#'
#' @description Estimate the variance for each block and then take the average.
#'
#' @example tests/testthat/examples/lm_variance.R
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
#' @rdname lm_variance
#' @export
lm_variance <- function(data, block_size = NULL, p = NULL) {
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
