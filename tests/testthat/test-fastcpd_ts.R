testthat::test_that(
  "invalid family provided for time series data", {
    testthat::expect_error(
      fastcpd.ts(
        data = seq_len(10),
        family = "at",
        order = 1
      ),
      r"[
The family should be one of "ar", "var" or "arima",
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
  "ar(1) model using `fastcpd_ts`", {
    set.seed(1)
    n <- 1000
    p <- 1
    x <- rep(0, n + 1)
    for (i in 1:600) {
      x[i + 1] <- 0.6 * x[i] + rnorm(1)
    }
    for (i in 601:1000) {
      x[i + 1] <- 0.3 * x[i] + rnorm(1)
    }

    result <- fastcpd_ts(x, "ar", 1)

    testthat::expect_equal(result@cp_set, 609)
  }
)

testthat::test_that(
  "ar(1) model using custom cost function", {
    set.seed(1)
    n <- 1000
    p <- 1
    time_series <- rep(0, n + 1)
    for (i in 1:600) {
      time_series[i + 1] <- 0.6 * time_series[i] + rnorm(1)
    }
    for (i in 601:1000) {
      time_series[i + 1] <- 0.3 * time_series[i] + rnorm(1)
    }
    ar1_loss <- function(data) {
      n <- nrow(data)
      optim_result <- optim(
        par = 0,
        fn = function(data, theta) {
          n <- nrow(data)
          unconditional_sum_of_squares <- (1 - theta^2) * data[1]^2 + sum(
            (data[2:n] - theta * data[1:(n - 1)])^2
          )
          log(unconditional_sum_of_squares) - log(1 - theta^2) / n
        },
        method = "Brent",
        data = data,
        lower = -0.999,
        upper = 0.999
      )
      n / 2 * (optim_result$value - log(n) + log(2 * pi) + 1)
    }
    result <- fastcpd(
      formula = ~ . - 1,
      data = data.frame(x = time_series[-1]),
      beta = (2 * p + 1) * log(n) / 2,
      p = 2 * p,
      cost = ar1_loss
    )

    testthat::expect_equal(result@cp_set, 614)
  }
)

testthat::test_that(
  "ma(4) model with `forecast::Arima` calculation of ML", {
    testthat::skip_if_not_installed("forecast")
    set.seed(1)
    n <- 400
    p <- 4
    time_series <- rep(0, n)

    lambda1 <- c(1, 3 / 2, 2i, -2i)
    param1 <- 1
    for (lam in lambda1) {
      param1 <- convolve(param1, c(1, -lam), type = "open")
    }
    ma_coef1 <- Re(param1[-1] / param1[1])

    testthat::expect_equal(ma_coef1, c(-5 / 3, 11 / 12, -5 / 12, 1 / 6))

    time_series[1:200] <-
      arima.sim(n = 200, model = list(ma = ma_coef1, sd = 0.1))

    lambda2 <- c(3, 1, -3i / 2, 3i / 2)
    param2 <- 1
    for (lam in lambda2) {
      param2 <- convolve(param2, c(1, -lam), type = "open")
    }
    ma_coef2 <- Re(param2[-1] / param2[1])

    testthat::expect_equal(ma_coef2, c(-4 / 3, 7 / 9, -16 / 27, 4 / 27))

    time_series[201:400] <-
      arima.sim(n = 200, model = list(ma = ma_coef2, sd = 0.3))

    result_ts <- fastcpd.ts(
      data = data.frame(x = time_series),
      family = "arima",
      order = c(0, 0, 4)
    )

    testthat::expect_equal(result_ts@cp_set, 200)

    result <- fastcpd(
      formula = ~ . - 1,
      data = data.frame(x = time_series),
      beta = (2 * p + 1) * log(n) / 2,
      order = c(0, 0, 4),
      family = "arima"
    )

    testthat::expect_equal(result@cp_set, 200)
  }
)

testthat::test_that(
  "example ar(3) model with innovation standard deviation 3", {
    set.seed(1)
    n <- 1000
    p <- 1
    x <- rep(0, n + 3)
    for (i in 1:600) {
      x[i + 3] <- 0.6 * x[i + 2] - 0.2 * x[i + 1] + 0.1 * x[i] + rnorm(1, 0, 3)
    }
    for (i in 601:1000) {
      x[i + 1] <- 0.3 * x[i + 2] + 0.4 * x[i + 1] + 0.2 * x[i] + rnorm(1, 0, 3)
    }

    result_ts <- fastcpd.ts(x, "ar", 3)

    testthat::expect_equal(result_ts@cp_set, 615)

    result <- fastcpd(
      formula = ~ . - 1,
      data = data.frame(x = x),
      order = 3,
      family = "ar"
    )

    testthat::expect_equal(result@cp_set, 615)
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

    result_ts <- fastcpd.ts(time_series[-1], "arima", c(1, 0, 0))

    testthat::expect_equal(result_ts@cp_set, 194)

    result <- fastcpd(
      formula = ~ . - 1,
      data = data.frame(x = time_series[-1]),
      beta = (2 * p + 1) * log(n) / 2,
      p = p,
      cost = ar1_loss
    )

    testthat::expect_equal(result@cp_set, 178)
  }
)

# testthat::test_that(
#   "garch model using rugarch package", {
#     set.seed(1)
#     new_env <- new.env()
#     load(
#       file = system.file("data", "dmbp.rda", package = "rugarch"),
#       envir = new_env
#     )
#     dmbp <- new_env[["dmbp"]]
#     result <- fastcpd(
#       formula = ~ x - 1,
#       data = data.frame(x = dmbp[, 1]),
#       cost = function(data) {
#         if (length(data) < 3) {
#           return(0)
#         }
#         spec <- rugarch::ugarchspec(
#           variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
#           mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
#           distribution.model = "norm"
#         )
#         fit <- rugarch::ugarchfit(spec = spec, data = data)
#         -fit@fit$LLH
#       }
#     )

#     testthat::expect_equal(result@cp_set, 100)
#   }
# )
