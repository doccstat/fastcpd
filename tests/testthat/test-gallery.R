# Please make sure this file is in sync with "vignettes/gallery.Rmd"

testthat::test_that(
  "linear regression with multi-dimensional responses", {
    testthat::skip_if_not_installed("mvtnorm")
    set.seed(1)
    n <- 300
    p <- 3
    y_count <- 2
    x <- mvtnorm::rmvnorm(n, rep(0, p), diag(p))
    theta_0 <- array(NA, dim = c(3, y_count, 3))
    theta_0[, , 1] <- cbind(c(1, 1.2, -1), c(-1, 0, 0.5))
    theta_0[, , 2] <- cbind(c(-1, 0, 0.5), c(0.5, -0.3, 0.2))
    theta_0[, , 3] <- cbind(c(0.5, -0.3, 0.2), c(1, 1.2, -1))
    y <- rbind(
      x[1:100, ] %*% theta_0[, , 1],
      x[101:200, ] %*% theta_0[, , 2],
      x[201:n, ] %*% theta_0[, , 3]
    ) + matrix(rnorm(n * y_count), ncol = y_count)
    multi_response_linear_loss <- function(data) {
      x <- data[, (ncol(data) - p + 1):ncol(data)]
      y <- data[, 1:(ncol(data) - p)]

      if (nrow(data) <= p) {
        x_t_x <- diag(p)
      } else {
        x_t_x <- crossprod(x)
      }

      norm(y - x %*% solve(x_t_x, t(x)) %*% y, type = "F")^2 / 2
    }
    result <- fastcpd(
      formula = y ~ x - 1,
      data = data.frame(y = y, x = x),
      beta = (2 * p + 1) * log(n) / 2,
      cost = multi_response_linear_loss,
      cp_only = TRUE
    )

    testthat::expect_equal(result@cp_set, c(102, 195))
  }
)

testthat::test_that(
  "custom cost function varince change in mean or var change data", {
    testthat::skip_if_not_installed("mvtnorm")
    set.seed(1)
    p <- 3
    data <- rbind.data.frame(
      mvtnorm::rmvnorm(300, mean = rep(0, p), sigma = diag(1, p)),
      mvtnorm::rmvnorm(400, mean = rep(10, p), sigma = diag(1, p)),
      mvtnorm::rmvnorm(300, mean = rep(0, p), sigma = diag(50, p)),
      mvtnorm::rmvnorm(300, mean = rep(0, p), sigma = diag(1, p)),
      mvtnorm::rmvnorm(400, mean = rep(10, p), sigma = diag(1, p)),
      mvtnorm::rmvnorm(300, mean = rep(10, p), sigma = diag(50, p))
    )
    data_all_mean <- colMeans(data)
    var_loss <- function(data) {
      n <- nrow(data)
      p <- ncol(data)
      if (n < p) {
        data_cov <- diag(p)
      } else {
        data_cov <- crossprod(sweep(data, 2, data_all_mean)) / (n - 1)
      }
      n / 2 * (log(det(data_cov)) + p * log(2 * pi) + p * (n - 1) / n)
    }
    var_loss_result <- fastcpd(
      formula = ~ . - 1,
      data = data,
      beta = (p^2 + 1) * log(nrow(data)) / 2,
      trim = 0.1,
      p = p^2,
      cost = var_loss
    )

    testthat::expect_equal(var_loss_result@cp_set, c(700, 1000, 1700))
  }
)

testthat::test_that(
  "example custom cost function huber regression", {
    testthat::skip_if_not_installed("mvtnorm")
    set.seed(1)
    n <- 400 + 300 + 500
    p <- 5
    x <- mvtnorm::rmvnorm(n, mean = rep(0, p), sigma = diag(p))
    theta <- rbind(
      mvtnorm::rmvnorm(1, mean = rep(0, p - 3), sigma = diag(p - 3)),
      mvtnorm::rmvnorm(1, mean = rep(5, p - 3), sigma = diag(p - 3)),
      mvtnorm::rmvnorm(1, mean = rep(9, p - 3), sigma = diag(p - 3))
    )
    theta <- cbind(theta, matrix(0, 3, 3))
    theta <- theta[rep(seq_len(3), c(400, 300, 500)), ]
    y_true <- rowSums(x * theta)
    factor <- c(
      2 * stats::rbinom(400, size = 1, prob = 0.95) - 1,
      2 * stats::rbinom(300, size = 1, prob = 0.95) - 1,
      2 * stats::rbinom(500, size = 1, prob = 0.95) - 1
    )
    y <- factor * y_true + stats::rnorm(n)
    data <- cbind.data.frame(y, x)
    huber_threshold <- 1
    huber_loss <- function(data, theta) {
      residual <- data[, 1] - data[, -1, drop = FALSE] %*% theta
      indicator <- abs(residual) <= huber_threshold
      sum(
        residual^2 / 2 * indicator +
          huber_threshold * (
            abs(residual) - huber_threshold / 2
          ) * (1 - indicator)
      )
    }
    huber_loss_gradient <- function(data, theta) {
      residual <- c(data[nrow(data), 1] - data[nrow(data), -1] %*% theta)
      if (abs(residual) <= huber_threshold) {
        -residual * data[nrow(data), -1]
      } else {
        -huber_threshold * sign(residual) * data[nrow(data), -1]
      }
    }
    huber_loss_hessian <- function(data, theta) {
      residual <- c(data[nrow(data), 1] - data[nrow(data), -1] %*% theta)
      if (abs(residual) <= huber_threshold) {
        outer(data[nrow(data), -1], data[nrow(data), -1])
      } else {
        0.01 * diag(length(theta))
      }
    }
    huber_regression_result <- fastcpd(
      formula = y ~ . - 1,
      data = data,
      beta = (p + 1) * log(n) / 2,
      cost = huber_loss,
      cost_gradient = huber_loss_gradient,
      cost_hessian = huber_loss_hessian
    )

    testthat::expect_equal(huber_regression_result@cp_set, c(401, 726))
  }
)

testthat::test_that(
  "huber regression with 0.1 vanilla", {
    testthat::skip_if_not_installed("mvtnorm")
    set.seed(1)
    n <- 400 + 300 + 500
    p <- 5
    x <- mvtnorm::rmvnorm(n, mean = rep(0, p), sigma = diag(p))
    theta <- rbind(
      mvtnorm::rmvnorm(1, mean = rep(0, p - 3), sigma = diag(p - 3)),
      mvtnorm::rmvnorm(1, mean = rep(5, p - 3), sigma = diag(p - 3)),
      mvtnorm::rmvnorm(1, mean = rep(9, p - 3), sigma = diag(p - 3))
    )
    theta <- cbind(theta, matrix(0, 3, 3))
    theta <- theta[rep(seq_len(3), c(400, 300, 500)), ]
    y_true <- rowSums(x * theta)
    factor <- c(
      2 * stats::rbinom(400, size = 1, prob = 0.95) - 1,
      2 * stats::rbinom(300, size = 1, prob = 0.95) - 1,
      2 * stats::rbinom(500, size = 1, prob = 0.95) - 1
    )
    y <- factor * y_true + stats::rnorm(n)
    data <- cbind.data.frame(y, x)
    huber_threshold <- 1
    huber_loss <- function(data, theta) {
      residual <- data[, 1] - data[, -1, drop = FALSE] %*% theta
      indicator <- abs(residual) <= huber_threshold
      sum(
        residual^2 / 2 * indicator +
          huber_threshold * (
            abs(residual) - huber_threshold / 2
          ) * (1 - indicator)
      )
    }
    huber_loss_gradient <- function(data, theta) {
      residual <- c(data[nrow(data), 1] - data[nrow(data), -1] %*% theta)
      if (abs(residual) <= huber_threshold) {
        -residual * data[nrow(data), -1]
      } else {
        -huber_threshold * sign(residual) * data[nrow(data), -1]
      }
    }
    huber_loss_hessian <- function(data, theta) {
      residual <- c(data[nrow(data), 1] - data[nrow(data), -1] %*% theta)
      if (abs(residual) <= huber_threshold) {
        outer(data[nrow(data), -1], data[nrow(data), -1])
      } else {
        0.01 * diag(length(theta))
      }
    }
    huber_regression_result <- fastcpd(
      formula = y ~ . - 1,
      data = data,
      beta = (p + 1) * log(n) / 2,
      cost = huber_loss,
      cost_gradient = huber_loss_gradient,
      cost_hessian = huber_loss_hessian,
      vanilla_percentage = 0.1
    )

    testthat::expect_equal(huber_regression_result@cp_set, c(401, 726))
  }
)
