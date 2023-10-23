testthat::test_that(
  "invalid family provided", {
    testthat::expect_error(
      fastcpd(
        data = data.frame(y = 0, x = 0),
        family = "bin0mial"
      ),
      r"[
The family should be one of
"gaussian", "binomial", "poisson", "lasso", "ar", "var", "arima", "garch",
"custom" or `NULL`, while the provided family is bin0mial.]"
    )
  }
)

testthat::test_that(
  "custom family but no cost function provided", {
    testthat::expect_error(
      fastcpd(
        data = data.frame(y = 0, x = 0)
      ),
      "cost function must be specified for custom family"
    )
  }
)

testthat::test_that(
  "linear regression without change points", {
    result <- fastcpd(
      formula = y ~ . - 1,
      data = data.frame(y = seq_len(100), x = seq_len(100)),
      family = "gaussian"
    )

    testthat::expect_length(result@cp_set, 0)
  }
)

testthat::test_that(
  "invalid family provided for time series data", {
    testthat::expect_error(
      fastcpd.ts(
        data = seq_len(10),
        family = "at",
        order = 1
      ),
      r"[
The family should be one of "ar", "var", "arima" or "garch",
while the provided family is at.]"
    )
  }
)

testthat::test_that(
  "invalid data provided for ar(1) model", {
    testthat::expect_error(
      fastcpd_ts(
        data = matrix(NA, 1, 2),
        family = "ar",
        order = 1
      ),
      "The data should be a univariate time series."
    )
  }
)

testthat::test_that(
  "invalid order for ar family", {
    testthat::expect_error(
      fastcpd(
        formula = ~ x - 1,
        data = data.frame(x = 0),
        family = "ar",
      ),
      "Please specify a positive integer as the `order` of the AR model."
    )
  }
)

testthat::test_that(
  "var(2) model with p = 2 with innovation sd = 1", {
    set.seed(1)
    n <- 1000
    p <- 2
    theta_1 <- matrix(c(-0.3, 0.6, -0.5, 0.4, 0.2, 0.1, 0.1, -0.2), nrow = p)
    theta_2 <- matrix(c(-0.1, -0.2, 0.1, 0.05, 0.1, -0.2, -0.2, 0.1), nrow = p)
    x <- matrix(0, n + 2, p)
    for (i in 1:600) {
      x[i + 2, ] <- theta_1 %*% c(x[i + 1, ], x[i, ]) + rnorm(p, 0, 1)
    }
    for (i in 601:1000) {
      x[i + 2, ] <- theta_2 %*% c(x[i + 1, ], x[i, ]) + rnorm(p, 0, 1)
    }
    data <- matrix(NA, n, p * 3)
    data[, seq_len(p)] <- x[p + seq_len(n), ]
    for (p_i in seq_len(p)) {
      data[, p + (p * p_i - 1):(p * p_i)] <- x[p - p_i + seq_len(n), ]
    }
    multi_response_linear_loss <- function(data) {
      x <- data[, (ncol(data) - p + 1):ncol(data)]
      y <- data[, 1:(ncol(data) - p)]

      if (nrow(data) <= p + 1) {
        x_t_x <- diag(p)
      } else {
        x_t_x <- crossprod(x)
      }

      norm(y - x %*% solve(x_t_x, t(x)) %*% y, type = "F")^2 / 2
    }
    result_custom <- fastcpd(
      formula = y.1 + y.2 ~ x.1 + x.2 + x.3 + x.4 - 1,
      data = data.frame(y = data[, 1:2], x = data[, 3:6]),
      beta = (2 * 4 + 1) * log(n) / 2,
      p = 4,
      cost = multi_response_linear_loss
    )

    testthat::expect_equal(result_custom@cp_set, 612)

    result_ts <- fastcpd_ts(x, "var", 2)

    # TODO(doccstat): Fix the order and p issues.
    testthat::expect_equal(result_ts@cp_set, c(424, 609, 720))

    result_var <- fastcpd(
      formula = ~ . - 1,
      data = data.frame(x = x),
      beta = (2 * 4 + 1) * log(n) / 2,
      order = 2,
      family = "var"
    )

    testthat::expect_equal(result_var@cp_set, 609)
  }
)

testthat::test_that(
  "ar(1) model with `forecast::Arima` calculation of ML", {
    testthat::skip_if_not_installed("forecast")
    set.seed(1)
    n <- 400
    p <- 1
    time_series <- rep(0, n + 1)
    for (i in 1:200) {
      time_series[i + 1] <- 0.6 * time_series[i] + rnorm(1)
    }
    for (i in 201:400) {
      time_series[i + 1] <- 0.3 * time_series[i] + rnorm(1)
    }
    ar1_loss <- function(data) {
      tryCatch(
        expr = -forecast::Arima(
          c(data), order = c(1, 0, 0), include.mean = FALSE, method = "ML"
        )$loglik,
        error = function(e) 0
      )
    }

    result <- fastcpd(
      formula = ~ . - 1,
      data = data.frame(x = time_series[-1]),
      beta = (2 * p + 1) * log(n) / 2,
      p = p,
      cost = ar1_loss,
      trim = 0
    )

    testthat::expect_equal(result@cp_set, c(177, 179))

    result_ts <-
      fastcpd.ts(time_series[-1], "arima", c(1, 0, 0), include.mean = FALSE)
    testthat::expect_equal(result_ts@cp_set, c(82, 178))
  }
)

testthat::test_that(
  "example linear regression with noise variance not equal to 1", {
    set.seed(1)
    p <- 4
    n <- 300
    cp <- c(100, 200)
    x <- mvtnorm::rmvnorm(n, rep(0, p), diag(p))
    theta_0 <- rbind(c(1, 3.2, -1, 0), c(-1, -0.5, 2.5, -2), c(0.8, -0.3, 1, 1))
    y <- c(
      x[1:cp[1], ] %*% theta_0[1, ] + rnorm(cp[1], 0, sd = 3),
      x[(cp[1] + 1):cp[2], ] %*% theta_0[2, ] + rnorm(cp[2] - cp[1], 0, sd = 3),
      x[(cp[2] + 1):n, ] %*% theta_0[3, ] + rnorm(n - cp[2], 0, sd = 3)
    )

    # Compute the variance estimation for each block and then take the average.
    block_size <- 5
    variance_estimation <- rep(NA, n - block_size)
    for (i in 1:(n - block_size)) {
      block_index <- seq_len(block_size) + i - 1
      block_index_lagged <- seq_len(block_size) + i

      y_block <- y[block_index]
      x_block <- x[block_index, ]

      y_block_lagged <- y[block_index_lagged]
      x_block_lagged <- x[block_index_lagged, ]

      x_t_x <- crossprod(x_block)
      x_t_x_lagged <- crossprod(x_block_lagged)

      block_slope <- solve(crossprod(x_block), crossprod(x_block, y_block))
      block_lagged_slope <- solve(
        crossprod(x_block_lagged), crossprod(x_block_lagged, y_block_lagged)
      )
      x_t_x_inv <- solve(x_t_x)
      x_t_x_inv_lagged <- solve(x_t_x_lagged)
      inv_product <- x_t_x_inv %*% x_t_x_inv_lagged
      cross_term <-
        inv_product %*% crossprod(x_block[-1, ], x_block_lagged[-block_size, ])
      delta_numerator <- crossprod(block_slope - block_lagged_slope)
      delta_denominator <-
        sum(diag(x_t_x_inv + x_t_x_inv_lagged - 2 * cross_term))
      variance_estimation[i] <- delta_numerator / delta_denominator
    }

    testthat::expect_equal(mean(variance_estimation), 10.472932)

    result <- fastcpd(
      data = data.frame(y = y, x = x),
      beta = (p + 1) * log(n) / 2 * mean(variance_estimation),
      family = "gaussian",
    )

    testthat::expect_equal(result@cp_set, c(100, 201))
  }
)
