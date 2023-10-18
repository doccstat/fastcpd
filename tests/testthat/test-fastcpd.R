testthat::test_that(
  "invalid family provided", {
    testthat::expect_error(
      fastcpd(
        data = data.frame(y = 0, x = 0),
        family = "bin0mial"
      ),
      r"[
The family should be one of "gaussian", "binomial", "poisson", "lasso", "ar",
"var", "custom" or `NULL` while the provided family is bin0mial.]"
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
  "invalid order for ar family", {
    testthat::expect_error(
      fastcpd(
        formula = ~ x - 1,
        data = data.frame(x = 0),
        family = "ar",
      ),
      "Please specify the order of the AR model as the parameter `p`."
    )
  }
)

testthat::test_that(
  "example linear regression", {
    testthat::skip_if_not_installed("mvtnorm")
    set.seed(1)
    p <- 3
    x <- mvtnorm::rmvnorm(300, rep(0, p), diag(p))
    theta_0 <- rbind(c(1, 1.2, -1), c(-1, 0, 0.5), c(0.5, -0.3, 0.2))
    y <- c(
      x[1:100, ] %*% theta_0[1, ] + rnorm(100, 0, 1),
      x[101:200, ] %*% theta_0[2, ] + rnorm(100, 0, 1),
      x[201:300, ] %*% theta_0[3, ] + rnorm(100, 0, 1)
    )
    result <- fastcpd(
      formula = y ~ . - 1,
      data = data.frame(y = y, x = x),
      family = "gaussian"
    )

    testthat::expect_equal(result@cp_set, c(98, 202))
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

    result_internal_variance <- fastcpd(
      data = data.frame(y = y, x = x),
      family = "gaussian",
    )

    testthat::expect_equal(
      result_internal_variance@cp_set, c(100, 201)
    )
  }
)

testthat::test_that(
  "example linear regression with one-dimensional covariate", {
    testthat::skip_if_not_installed("mvtnorm")
    set.seed(1)
    p <- 1
    x <- mvtnorm::rmvnorm(300, rep(0, p), diag(p))
    theta_0 <- matrix(c(1, -1, 0.5))
    y <- c(
      x[1:100, ] * theta_0[1, ] + rnorm(100, 0, 1),
      x[101:200, ] * theta_0[2, ] + rnorm(100, 0, 1),
      x[201:300, ] * theta_0[3, ] + rnorm(100, 0, 1)
    )
    result <- fastcpd(
      formula = y ~ . - 1,
      data = data.frame(y = y, x = x),
      family = "gaussian"
    )

    testthat::expect_equal(result@cp_set, c(100, 194))
  }
)

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
  "example logistic regression", {
    set.seed(1)
    x <- matrix(rnorm(1500, 0, 1), ncol = 5)
    theta <- rbind(rnorm(5, 0, 1), rnorm(5, 2, 1))
    y <- c(
      rbinom(125, 1, 1 / (1 + exp(-x[1:125, ] %*% theta[1, ]))),
      rbinom(175, 1, 1 / (1 + exp(-x[126:300, ] %*% theta[2, ])))
    )
    result <- suppressWarnings(fastcpd(
      formula = y ~ . - 1,
      data = data.frame(y = y, x = x),
      family = "binomial"
    ))

    testthat::expect_equal(result@cp_set, 126)
  }
)

testthat::test_that(
  "logistic regression with warm start", {
    set.seed(1)
    x <- matrix(rnorm(1500, 0, 1), ncol = 5)
    theta <- rbind(rnorm(5, 0, 1), rnorm(5, 2, 1))
    y <- c(
      rbinom(125, 1, 1 / (1 + exp(-x[1:125, ] %*% theta[1, ]))),
      rbinom(175, 1, 1 / (1 + exp(-x[126:300, ] %*% theta[2, ])))
    )
    result <- suppressWarnings(fastcpd(
      formula = y ~ . - 1,
      data = data.frame(y = y, x = x),
      family = "binomial",
      vanilla_percentage = 1,
      warm_start = TRUE
    ))

    testthat::expect_equal(result@cp_set, 121)
  }
)

testthat::test_that(
  "example poisson regression", {
    testthat::skip_if_not_installed("mvtnorm")
    set.seed(1)
    p <- 3
    x <- mvtnorm::rmvnorm(1500, rep(0, p), diag(p))
    delta <- rnorm(p)
    theta_0 <- c(1, 1.2, -1)
    y <- c(
      rpois(300, exp(x[1:300, ] %*% theta_0)),
      rpois(400, exp(x[301:700, ] %*% (theta_0 + delta))),
      rpois(300, exp(x[701:1000, ] %*% theta_0)),
      rpois(100, exp(x[1001:1100, ] %*% (theta_0 - delta))),
      rpois(200, exp(x[1101:1300, ] %*% theta_0)),
      rpois(200, exp(x[1301:1500, ] %*% (theta_0 + delta)))
    )
    result <- fastcpd(
      formula = y ~ . - 1,
      data = data.frame(y = y, x = x),
      beta = (p + 1) * log(1500) / 2,
      k = function(x) 0,
      family = "poisson",
      epsilon = 1e-5
    )

    result_two_epochs <- fastcpd(
      formula = y ~ . - 1,
      data = data.frame(y = y, x = x),
      beta = (p + 1) * log(1500) / 2,
      k = function(x) 1,
      family = "poisson",
      epsilon = 1e-4
    )

    testthat::expect_equal(result@cp_set, c(329, 728, 1021, 1107, 1325))
    testthat::expect_equal(
      result_two_epochs@cp_set, c(328, 716, 1020, 1102, 1323)
    )
  }
)

testthat::test_that(
  "example penalized linear regression", {
    testthat::skip_if_not_installed("mvtnorm")
    set.seed(1)
    n <- 1500
    p_true <- 6
    p <- 50
    x <- mvtnorm::rmvnorm(1500, rep(0, p), diag(p))
    theta_0 <- rbind(
      runif(p_true, -5, -2),
      runif(p_true, -3, 3),
      runif(p_true, 2, 5),
      runif(p_true, -5, 5)
    )
    theta_0 <- cbind(theta_0, matrix(0, ncol = p - p_true, nrow = 4))
    y <- c(
      x[1:300, ] %*% theta_0[1, ] + rnorm(300, 0, 1),
      x[301:700, ] %*% theta_0[2, ] + rnorm(400, 0, 1),
      x[701:1000, ] %*% theta_0[3, ] + rnorm(300, 0, 1),
      x[1001:1500, ] %*% theta_0[4, ] + rnorm(500, 0, 1)
    )
    result <- fastcpd(
      formula = y ~ . - 1,
      data = data.frame(y = y, x = x),
      family = "lasso"
    )

    testthat::expect_equal(result@cp_set, c(300, 700, 1000))
    testthat::expect_equal(
      unname(sqrt(colSums((t(theta_0[, 1:6]) - result@thetas[1:6, ])^2))),
      c(0.5801021, 0.5439067, 0.4677495, 0.4898553),
      tolerance = 1e-7,
      ignore_attr = TRUE
    )
  }
)

testthat::test_that(
  "penalized linear regression with smaller sample size", {
    testthat::skip_if_not_installed("mvtnorm")
    set.seed(1)
    n <- 1500
    p_true <- 6
    p <- 50
    x <- mvtnorm::rmvnorm(480, rep(0, p), diag(p))
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
      x[321:480, ] %*% theta_0[4, ] + rnorm(160, 0, 1)
    )
    result <- fastcpd(
      formula = y ~ . - 1,
      data = data.frame(y = y, x = x),
      family = "lasso"
    )

    testthat::expect_equal(result@cp_set, c(79, 201, 325))

    result_multiple_epochs <- fastcpd(
      formula = y ~ . - 1,
      data = data.frame(y = y, x = x),
      family = "lasso",
      k = function(x) if (x < 20) 1 else 0
    )

    testthat::expect_equal(result_multiple_epochs@cp_set, c(80, 200, 320))

    result_vanilla_percentage <- fastcpd(
      formula = y ~ . - 1,
      data = data.frame(y = y, x = x),
      family = "lasso",
      vanilla_percentage = 0.2
    )

    testthat::expect_equal(result_vanilla_percentage@cp_set, c(79, 201, 321))
  }
)

testthat::test_that(
  "example custom logistic regression", {
    set.seed(1)
    p <- 5
    x <- matrix(rnorm(375 * p, 0, 1), ncol = p)
    theta <- rbind(rnorm(p, 0, 1), rnorm(p, 2, 1))
    y <- c(
      rbinom(200, 1, 1 / (1 + exp(-x[1:200, ] %*% theta[1, ]))),
      rbinom(175, 1, 1 / (1 + exp(-x[201:375, ] %*% theta[2, ])))
    )
    data <- data.frame(y = y, x = x)
    warning_messages <- testthat::capture_warnings(
      result_builtin <- fastcpd(
        formula = y ~ . - 1,
        data = data,
        family = "binomial"
      )
    )
    logistic_loss <- function(data, theta) {
      x <- data[, -1]
      y <- data[, 1]
      u <- x %*% theta
      nll <- -y * u + log(1 + exp(u))
      nll[u > 10] <- -y[u > 10] * u[u > 10] + u[u > 10]
      sum(nll)
    }
    logistic_loss_gradient <- function(data, theta) {
      x <- data[nrow(data), -1]
      y <- data[nrow(data), 1]
      c(-(y - 1 / (1 + exp(-x %*% theta)))) * x
    }
    logistic_loss_hessian <- function(data, theta) {
      x <- data[nrow(data), -1]
      prob <- 1 / (1 + exp(-x %*% theta))
      (x %o% x) * c((1 - prob) * prob)
    }
    result_custom <- fastcpd(
      formula = y ~ . - 1,
      data = data,
      epsilon = 1e-5,
      cost = logistic_loss,
      cost_gradient = logistic_loss_gradient,
      cost_hessian = logistic_loss_hessian
    )

    result_custom_two_epochs <- fastcpd(
      formula = y ~ . - 1,
      data = data,
      k = function(x) 1,
      epsilon = 1e-5,
      cost = logistic_loss,
      cost_gradient = logistic_loss_gradient,
      cost_hessian = logistic_loss_hessian
    )

    testthat::expect_equal(
      warning_messages,
      rep("fit_glm: fitted probabilities numerically 0 or 1 occurred", 3)
    )

    testthat::expect_equal(result_builtin@cp_set, 200)
    testthat::expect_equal(result_custom@cp_set, 201)
    testthat::expect_equal(result_custom_two_epochs@cp_set, 200)
  }
)

testthat::test_that(
  "custom one-dimensional logistic regression", {
    set.seed(1)
    p <- 1
    x <- matrix(rnorm(375 * p, 0, 1), ncol = p)
    theta <- rbind(rnorm(p, 0, 1), rnorm(p, 2, 1))
    y <- c(
      rbinom(200, 1, 1 / (1 + exp(-x[1:200, ] * theta[1, ]))),
      rbinom(175, 1, 1 / (1 + exp(-x[201:375, ] * theta[2, ])))
    )
    data <- data.frame(y = y, x = x)

    result_builtin <- fastcpd(
      formula = y ~ . - 1,
      data = data,
      family = "binomial"
    )
    logistic_loss <- function(data, theta) {
      x <- data[, -1]
      y <- data[, 1]
      u <- x * c(theta)
      nll <- -y * u + log(1 + exp(u))
      nll[u > 10] <- -y[u > 10] * u[u > 10] + u[u > 10]
      sum(nll)
    }
    logistic_loss_gradient <- function(data, theta) {
      x <- data[nrow(data), -1]
      y <- data[nrow(data), 1]
      c(-(y - 1 / (1 + exp(-x * theta)))) * x
    }
    logistic_loss_hessian <- function(data, theta) {
      x <- data[nrow(data), -1]
      prob <- 1 / (1 + exp(-x * theta))
      (x %o% x) * c((1 - prob) * prob)
    }
    result_custom <- fastcpd(
      formula = y ~ . - 1,
      data = data,
      epsilon = 1e-5,
      cost = logistic_loss,
      cost_gradient = logistic_loss_gradient,
      cost_hessian = logistic_loss_hessian
    )

    testthat::expect_equal(result_builtin@cp_set, 198)
    testthat::expect_equal(result_custom@cp_set, 200)
  }
)

testthat::test_that(
  "example custom cost function mean change", {
    testthat::skip_if_not_installed("mvtnorm")
    set.seed(1)
    p <- 1
    data <- rbind(
      mvtnorm::rmvnorm(300, mean = rep(0, p), sigma = diag(100, p)),
      mvtnorm::rmvnorm(400, mean = rep(50, p), sigma = diag(100, p)),
      mvtnorm::rmvnorm(300, mean = rep(2, p), sigma = diag(100, p))
    )
    segment_count_guess <- 10
    block_size <- max(floor(sqrt(nrow(data)) / (segment_count_guess + 1)), 2)
    block_count <- floor(nrow(data) / block_size)
    data_all_vars <- rep(0, block_count)
    for (block_index in seq_len(block_count)) {
      block_start <- (block_index - 1) * block_size + 1
      block_end <- if (block_index < block_count) {
        block_index * block_size
      } else {
        nrow(data)
      }
      data_all_vars[block_index] <- var(data[block_start:block_end, ])
    }
    data_all_var <- mean(data_all_vars)
    mean_loss <- function(data) {
      n <- nrow(data)
      n / 2 * (
        log(data_all_var) + log(2 * pi) +
          sum((data - colMeans(data))^2 / data_all_var) / n
      )
    }
    mean_loss_result <- fastcpd(
      formula = ~ . - 1,
      data = data.frame(data),
      beta = (p + 1) * log(nrow(data)) / 2,
      p = p,
      cost = mean_loss
    )

    testthat::expect_equal(mean_loss_result@cp_set, c(300, 700))
  }
)

testthat::test_that(
  "example custom cost function multivariate mean change", {
    testthat::skip_if_not_installed("mvtnorm")
    set.seed(1)
    p <- 3
    data <- rbind(
      mvtnorm::rmvnorm(300, mean = rep(0, p), sigma = diag(100, p)),
      mvtnorm::rmvnorm(400, mean = rep(50, p), sigma = diag(100, p)),
      mvtnorm::rmvnorm(300, mean = rep(2, p), sigma = diag(100, p))
    )
    segment_count_guess <- 5
    block_size <- max(floor(sqrt(nrow(data)) / (segment_count_guess + 1)), 2)
    block_count <- floor(nrow(data) / block_size)
    data_all_covs <- array(NA, dim = c(block_count, p, p))
    for (block_index in seq_len(block_count)) {
      block_start <- (block_index - 1) * block_size + 1
      block_end <- if (block_index < block_count) {
        block_index * block_size
      } else {
        nrow(data)
      }
      data_all_covs[block_index, , ] <- cov(data[block_start:block_end, ])
    }
    data_all_cov <- colMeans(data_all_covs)
    mean_loss <- function(data) {
      n <- nrow(data)
      demeaned_data <- sweep(data, 2, colMeans(data))
      n / 2 * (
        log(det(data_all_cov)) + p * log(2 * pi) +
          sum(diag(solve(data_all_cov, crossprod(demeaned_data)))) / n
      )
    }
    mean_loss_result <- fastcpd(
      formula = ~ . - 1,
      data = data.frame(data),
      beta = (p + 1) * log(nrow(data)) / 2,
      p = p,
      cost = mean_loss
    )

    testthat::expect_equal(mean_loss_result@cp_set, c(300, 700))
  }
)

testthat::test_that(
  "example custom cost function variance change", {
    testthat::skip_if_not_installed("mvtnorm")
    set.seed(1)
    p <- 1
    data <- rbind.data.frame(
      mvtnorm::rmvnorm(300, mean = rep(0, p), sigma = diag(1, p)),
      mvtnorm::rmvnorm(400, mean = rep(0, p), sigma = diag(50, p)),
      mvtnorm::rmvnorm(300, mean = rep(0, p), sigma = diag(2, p))
    )
    data_all_mean <- colMeans(data)
    var_loss <- function(data) {
      n <- nrow(data)
      data_cov <- crossprod(sweep(data, 2, data_all_mean)) / (n - 1)
      n / 2 * (log(data_cov) + log(2 * pi) + (n - 1) / n)
    }
    var_loss_result <- fastcpd(
      formula = ~ . - 1,
      data = data,
      beta = (p + 1) * log(nrow(data)) / 2,
      p = p,
      cost = var_loss
    )

    testthat::expect_equal(var_loss_result@cp_set, c(300, 699))
  }
)

testthat::test_that(
  "example custom cost function multivariate variance change", {
    testthat::skip_if_not_installed("mvtnorm")
    set.seed(1)
    p <- 3
    data <- rbind.data.frame(
      mvtnorm::rmvnorm(
        300, rep(0, p), crossprod(matrix(runif(p^2) * 2 - 1, p))
      ),
      mvtnorm::rmvnorm(
        400, rep(0, p), crossprod(matrix(runif(p^2) * 2 - 1, p))
      ),
      mvtnorm::rmvnorm(
        300, rep(0, p), crossprod(matrix(runif(p^2) * 2 - 1, p))
      )
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

    testthat::expect_equal(var_loss_result@cp_set, c(300, 700))
  }
)

testthat::test_that(
  "example custom cost function mean or variance change", {
    testthat::skip_if_not_installed("mvtnorm")
    set.seed(1)
    p <- 1
    data <- rbind.data.frame(
      mvtnorm::rmvnorm(300, mean = rep(0, p), sigma = diag(1, p)),
      mvtnorm::rmvnorm(400, mean = rep(10, p), sigma = diag(1, p)),
      mvtnorm::rmvnorm(300, mean = rep(0, p), sigma = diag(50, p)),
      mvtnorm::rmvnorm(300, mean = rep(0, p), sigma = diag(1, p)),
      mvtnorm::rmvnorm(400, mean = rep(10, p), sigma = diag(1, p)),
      mvtnorm::rmvnorm(300, mean = rep(10, p), sigma = diag(50, p))
    )
    meanvar_loss <- function(data) {
      n <- nrow(data)
      data_cov <- 1
      if (n > 1) {
        data_cov <- var(data)
      }
      n / 2 * (log(data_cov) + log(2 * pi) + (n - 1) / n)
    }
    meanvar_loss_result <- fastcpd(
      formula = ~ . - 1,
      data = data,
      beta = (p^2 + p + 1) * log(nrow(data)) / 2,
      p = p^2 + p,
      cost = meanvar_loss
    )

    testthat::expect_equal(
      meanvar_loss_result@cp_set, c(300, 700, 1000, 1300, 1700)
    )
  }
)

testthat::test_that(
  "example custom cost function multivariate mean or variance change", {
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
    meanvar_loss <- function(data) {
      n <- nrow(data)
      p <- ncol(data)
      if (n <= p) {
        data_cov <- diag(p)
      } else {
        data_cov <- cov(data)
      }
      n / 2 * (log(det(data_cov)) + p * log(2 * pi) + p * (n - 1) / n)
    }
    meanvar_loss_result <- fastcpd(
      formula = ~ . - 1,
      data = data,
      beta = (p^2 + p + 1) * log(nrow(data)) / 2,
      trim = 0.01,
      p = p^2 + p,
      cost = meanvar_loss
    )

    testthat::expect_equal(
      meanvar_loss_result@cp_set, c(300, 700, 1000, 1300, 1700)
    )
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

testthat::test_that(
  "example ar(1) model", {
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

    result <- fastcpd(
      formula = ~ . - 1,
      data = data.frame(x = x),
      p = 1,
      family = "ar"
    )

    testthat::expect_equal(result@cp_set, 609)
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

    result <- fastcpd(
      formula = ~ . - 1,
      data = data.frame(x = x),
      p = 3,
      family = "ar"
    )

    testthat::expect_equal(result@cp_set, 615)
  }
)

testthat::test_that(
  "example var(2) model with p = 2 with innovation sd = 1", {
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
    result <- fastcpd(
      formula = y.1 + y.2 ~ x.1 + x.2 + x.3 + x.4 - 1,
      data = data.frame(y = data[, 1:2], x = data[, 3:6]),
      beta = (2 * 4 + 1) * log(n) / 2,
      p = 4,
      cost = multi_response_linear_loss
    )

    testthat::expect_equal(result@cp_set, 612)

    result_var <- fastcpd(
      formula = ~ . - 1,
      data = data.frame(x = x),
      beta = (2 * 4 + 1) * log(n) / 2,
      p = 2,
      family = "var"
    )

    testthat::expect_equal(result_var@cp_set, 611)
  }
)

testthat::test_that(
  "ar(1) model", {
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
        error = function(e) {
          0
        }
      )
    }
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

testthat::test_that(
  "var(1) model with p = 2", {
    set.seed(1)
    n <- 1000
    p <- 2
    theta_0 <- array(NA, dim = c(p, p, 2))
    theta_0[, , 1] <- cbind(c(-0.3, 0.6), c(-0.4, 0.5))
    theta_0[, , 2] <- cbind(c(-0.1, -0.2), c(0.1, 0.05))
    x <- matrix(0, n + 1, p)
    for (i in 1:600) {
      x[i + 1, ] <- theta_0[, , 1] %*% x[i, ] + rnorm(p)
    }
    for (i in 601:1000) {
      x[i + 1, ] <- theta_0[, , 2] %*% x[i, ] + rnorm(p)
    }
    y <- x[3:1001, ]
    x <- x[2:1000, ]
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
      cost = multi_response_linear_loss
    )

    testthat::expect_equal(result@cp_set, 593)
  }
)
