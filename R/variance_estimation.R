#' @title Variance estimation for ARMA model with change points
#'
#' @description Estimate the variance for each block and then take the average.
#'
#' @example tests/testthat/examples/variance_arma.R
#'
#' @md
#'
#' @param data A one-column matrix or a vector.
#' @param p The order of the autoregressive part.
#' @param q The order of the moving average part.
#'
#' @return A numeric value representing the variance.
#'
#' @rdname variance_arma
#' @export
variance_arma <- function(data, p, q) {
  p <- min(3 * p, 15)
  y <- data[p + seq_len(length(data) - p)]
  x <- matrix(NA, length(data) - p, p)
  for (p_i in seq_len(p)) {
    x[, p_i] <- data[(p - p_i) + seq_len(length(data) - p)]
  }
  variance.lm(cbind(y, x), outlier_iqr = 3)
}

#' @rdname variance_arma
#' @export
variance.arma <- variance_arma  # nolint: Conventional R function style

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
#' @param d The dimension of the response variable.
#' @param block_size The size of the blocks to use for variance estimation.
#'
#' @return A numeric value representing the variance.
#'
#' @rdname variance_lm
#' @export
variance_lm <- function(
  data,
  d = 1,
  block_size = ncol(data) - d + 1,
  outlier_iqr = Inf
) {
  data <- as.matrix(data)
  n <- nrow(data)
  estimators <- array(NA, c(n - block_size, d, d))
  for (i in seq_len(n - block_size)) {
    block_index <- seq_len(block_size) + i - 1
    block_index_lagged <- seq_len(block_size) + i

    y_block <- data[block_index, seq_len(d), drop = FALSE]
    x_block <- data[block_index, -seq_len(d), drop = FALSE]

    y_block_lagged <- data[block_index_lagged, seq_len(d), drop = FALSE]
    x_block_lagged <- data[block_index_lagged, -seq_len(d), drop = FALSE]

    x_t_x <- crossprod(x_block)
    x_t_x_lagged <- crossprod(x_block_lagged)

    tryCatch(
      expr = {
        block_slope <- solve(crossprod(x_block), crossprod(x_block, y_block))
        block_lagged_slope <- solve(
          crossprod(x_block_lagged),
          crossprod(x_block_lagged, y_block_lagged)
        )
        x_t_x_inv <- solve(x_t_x)
        x_t_x_inv_lagged <- solve(x_t_x_lagged)
        cross_term_x <- crossprod(
          x_block[-1, , drop = FALSE],
          x_block_lagged[-block_size, , drop = FALSE]
        )
        cross_term <- x_t_x_inv %*% x_t_x_inv_lagged %*% cross_term_x
        delta_numerator <- crossprod(block_slope - block_lagged_slope)
        delta_denominator <- matrix(0, d, d)
        for (j in seq_len(d)) {
          for (k in seq_len(d)) {
            if (j != k) {
              delta_denominator[j, k] <- delta_denominator[j, k] + crossprod(
                block_slope[, j] - block_lagged_slope[, k]
              )
            }
          }
        }
        delta_denominator <- delta_denominator + sum(
          diag(x_t_x_inv + x_t_x_inv_lagged - 2 * cross_term)
        )
        estimators[i, , ] <- delta_numerator / delta_denominator
      },
      error = function(e) {
        estimators[i, , ] <- matrix(NA, d, d)
        message("Variance estimation failed for block ", i, ".")
      }
    )
  }
  if (d == 1) {
    estimators <- stats::na.exclude(c(estimators))
    outlier_threshold <-
      stats::quantile(estimators, 0.75) + outlier_iqr * stats::IQR(estimators)
    mean(estimators[estimators < outlier_threshold], na.rm = TRUE)
  } else {
    colMeans(estimators, na.rm = TRUE)
  }
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
  estimators <- array(NA, c(n - 1, p, p))
  for (i in seq_len(n - 1)) {
    estimators[i, , ] <- crossprod(
      data[i + 1, , drop = FALSE] - data[i, , drop = FALSE]
    )
  }
  colMeans(estimators, na.rm = TRUE) / 2
}

#' @rdname variance_mean
#' @export
variance.mean <- variance_mean  # nolint: Conventional R function style
