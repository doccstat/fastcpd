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
