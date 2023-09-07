testthat::test_that("linear regression with one-dimensional covariate", {
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
})

testthat::test_that("random linear regression", {
  result <- fastcpd(
    formula = y ~ . - 1,
    data = data.frame(y = seq_len(100), x = seq_len(100)),
    family = "gaussian"
  )

  testthat::expect_length(result@cp_set, 0)
})

testthat::test_that("logistic regression", {
  # This is the same example with `fastcpd` documentation. Please keep it in
  # sync if the documentation ever changes.
  set.seed(1)

  kChangePointLocation <- 125  # nolint: Google Style Guide
  kNumberOfDataPoints <- 300  # nolint: Google Style Guide
  kDimension <- 5  # nolint: Google Style Guide

  # There are 300 five-dimensional data points.
  x <- matrix(rnorm(kNumberOfDataPoints * kDimension, 0, 1), ncol = kDimension)

  # Randomly generate coefficients with different means.
  theta <- rbind(rnorm(kDimension, 0, 1), rnorm(kDimension, 2, 1))

  # Randomly generate response variables based on the segmented data and
  # corresponding coefficients
  y <- c(
    rbinom(
      kChangePointLocation,
      1,
      1 / (1 + exp(-x[1:kChangePointLocation, ] %*% theta[1, ]))
    ),
    rbinom(
      kNumberOfDataPoints - kChangePointLocation,
      1,
      1 / (1 + exp(
        -x[(kChangePointLocation + 1):kNumberOfDataPoints, ] %*% theta[2, ]
      ))
    )
  )

  change_points_binomial_fastcpd <- fastcpd(
    formula = y ~ . - 1,
    data = data.frame(y = y, x = x),
    family = "binomial",
    segment_count = 5
  )@cp_set

  testthat::expect_equal(change_points_binomial_fastcpd, kChangePointLocation)
})

testthat::test_that("linear regression", {
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
})

testthat::test_that("logistic regression", {
  # This is the same example with `fastcpd` documentation. Please keep it in
  # sync if the documentation ever changes.
  set.seed(1)

  kChangePointLocation <- 125  # nolint: Google Style Guide
  kNumberOfDataPoints <- 300  # nolint: Google Style Guide
  kDimension <- 5  # nolint: Google Style Guide

  # There are 300 five-dimensional data points.
  x <- matrix(rnorm(kNumberOfDataPoints * kDimension, 0, 1), ncol = kDimension)

  # Randomly generate coefficients with different means.
  theta <- rbind(rnorm(kDimension, 0, 1), rnorm(kDimension, 2, 1))

  # Randomly generate response variables based on the segmented data and
  # corresponding coefficients
  y <- c(
    rbinom(
      kChangePointLocation,
      1,
      1 / (1 + exp(-x[1:kChangePointLocation, ] %*% theta[1, ]))
    ),
    rbinom(kNumberOfDataPoints - kChangePointLocation, 1, 1 / (1 + exp(
      -x[(kChangePointLocation + 1):kNumberOfDataPoints, ] %*% theta[2, ]
    )))
  )

  change_points_binomial_fastcpd <- fastcpd(
    formula = y ~ . - 1,
    data = data.frame(y = y, x = x),
    family = "binomial",
    segment_count = 5
  )@cp_set

  testthat::expect_equal(change_points_binomial_fastcpd, kChangePointLocation)
})

testthat::test_that("poisson regression", {
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
})

testthat::test_that("penalized linear regression", {
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

  testthat::expect_equal(result@cp_set, c(300, 701, 1000))
})

testthat::test_that("custom logistic regression", {
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

  testthat::expect_equal(
    warning_messages,
    rep("fit_glm: fitted probabilities numerically 0 or 1 occurred", 3)
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

  testthat::expect_equal(result_builtin@cp_set, 200)
  testthat::expect_equal(result_custom@cp_set, 201)
  testthat::expect_equal(result_custom_two_epochs@cp_set, 200)
})

testthat::test_that("custom one-dimensional logistic regression", {
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
})

testthat::test_that("mean change", {
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
    (norm(data, type = "F")^2 - colSums(data)^2 / n) / 2 / data_all_var +
      n / 2 * (log(data_all_var) + log(2 * pi))
  }
  mean_loss_result <- fastcpd(
    formula = ~ . - 1,
    data = data.frame(data),
    beta = (p + 1) * log(nrow(data)) / 2,
    p = p,
    cost = mean_loss
  )

  testthat::expect_equal(mean_loss_result@cp_set, c(300, 700))
})

testthat::test_that("variance change", {
  set.seed(1)
  p <- 1
  data <- rbind.data.frame(
    mvtnorm::rmvnorm(300, mean = rep(0, p), sigma = diag(1, p)),
    mvtnorm::rmvnorm(400, mean = rep(0, p), sigma = diag(50, p)),
    mvtnorm::rmvnorm(300, mean = rep(0, p), sigma = diag(2, p))
  )
  data_all_mu <- colMeans(data)
  var_loss <- function(data) {
    demeaned_data_norm <- norm(sweep(data, 2, data_all_mu), type = "F")
    nrow(data) * (1 + log(2 * pi) + log(demeaned_data_norm^2 / nrow(data))) / 2
  }
  var_loss_result <- fastcpd(
    formula = ~ . - 1,
    data = data,
    beta = (p + 1) * log(nrow(data)) / 2,
    p = p,
    cost = var_loss
  )

  testthat::expect_equal(var_loss_result@cp_set, c(300, 699))
})

testthat::test_that("mean / variance change", {
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
    loss_part <- (colSums(data^2) - colSums(data)^2 / nrow(data)) / nrow(data)
    nrow(data) * (1 + log(2 * pi) + log(loss_part)) / 2
  }
  meanvar_loss_result <- fastcpd(
    formula = ~ . - 1,
    data = data,
    beta = (2 * p + 1) * log(nrow(data)) / 2,
    p = 2 * p,
    cost = meanvar_loss
  )

  testthat::expect_equal(
    meanvar_loss_result@cp_set, c(300, 700, 1000, 1300, 1700)
  )
})

testthat::test_that("huber regression", {
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
})

testthat::test_that("huber regression with 0.1 vanilla", {
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
})

testthat::test_that("confidence interval experiment", {
  testthat::skip("This test is intended to be run manually.")
  set.seed(1)
  kDimension <- 1  # nolint: Google Style Guide
  change_point_locations <- NULL
  for (experiment_id in seq_len(20)) {
    data <- rbind(
      mvtnorm::rmvnorm(
        300, mean = rep(0, kDimension), sigma = diag(100, kDimension)
      ),
      mvtnorm::rmvnorm(
        400, mean = rep(50, kDimension), sigma = diag(100, kDimension)
      ),
      mvtnorm::rmvnorm(
        300, mean = rep(2, kDimension), sigma = diag(100, kDimension)
      )
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
      (norm(data, type = "F")^2 - colSums(data)^2 / n) / 2 / data_all_var +
        n / 2 * (log(data_all_var) + log(2 * pi))
    }
    mean_loss_result <- fastcpd(
      formula = ~ . - 1,
      data = data.frame(data),
      beta = (kDimension + 1) * log(nrow(data)) / 2,
      p = kDimension,
      cost = mean_loss
    )
    change_point_locations <- c(change_point_locations, mean_loss_result@cp_set)
  }

  change_locations_cookie_bucket <- NULL
  cookie_bucket_id_list <- sample.int(n = 20, size = 1000, replace = TRUE)
  all_data <- data
  for (cookie_bucket_id in seq_len(20)) {
    data <- all_data[cookie_bucket_id_list != cookie_bucket_id, , drop = FALSE]
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
      (norm(data, type = "F")^2 - colSums(data)^2 / n) / 2 / data_all_var +
        n / 2 * (log(data_all_var) + log(2 * pi))
    }
    mean_loss_result <- fastcpd(
      formula = ~ . - 1,
      data = data.frame(data),
      beta = (kDimension + 1) * log(nrow(data)) / 2,
      p = kDimension,
      cost = mean_loss
    )

    for (cp in mean_loss_result@cp_set) {
      ordinal_mapped_cp <- which(
        cumsum(cookie_bucket_id_list != cookie_bucket_id) == cp
      )[1]
      change_locations_cookie_bucket <-
        c(change_locations_cookie_bucket, ordinal_mapped_cp)
    }
  }

  table_change_point_locations <- table(change_point_locations)
  testthat::expect_equal(
    rownames(table_change_point_locations), c("269", "299", "300", "700")
  )
  testthat::expect_equal(
    unname(table_change_point_locations), c(1, 1, 19, 20), ignore_attr = TRUE
  )

  table_cp_cookie_bucket <- table(change_locations_cookie_bucket)
  testthat::expect_equal(
    rownames(table_cp_cookie_bucket), c("299", "300", "697", "700")
  )
  testthat::expect_equal(
    unname(table_cp_cookie_bucket), c(1, 19, 1, 19), ignore_attr = TRUE
  )
})

testthat::test_that("confidence interval experiment with one change point", {
  testthat::skip("This test is intended to be run manually.")
  set.seed(1)
  kDimension <- 1  # nolint: Google Style Guide
  change_point_locations <- list()
  change_locations_cookie_bucket <- rep(list(rep(list(NULL), 20)), 500)
  containing_change_point <- rep(FALSE, 500)
  for (experiment_id in seq_len(500)) {
    data <- rbind(
      mvtnorm::rmvnorm(
        121, mean = rep(0, kDimension), sigma = diag(100, kDimension)
      ),
      mvtnorm::rmvnorm(
        121, mean = rep(50, kDimension), sigma = diag(100, kDimension)
      )
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
      (norm(data, type = "F")^2 - colSums(data)^2 / n) / 2 / data_all_var +
        n / 2 * (log(data_all_var) + log(2 * pi))
    }
    mean_loss_result <- fastcpd(
      formula = ~ . - 1,
      data = data.frame(data),
      beta = (kDimension + 1) * log(nrow(data)) / 2,
      p = kDimension,
      cost = mean_loss
    )
    change_point_locations[[experiment_id]] <- mean_loss_result@cp_set

    cookie_bucket_id_list <-
      sample.int(n = 20, size = 121 + 121, replace = TRUE)
    all_data <- data
    for (cookie_bucket_id in seq_len(20)) {
      data <-
        all_data[cookie_bucket_id_list != cookie_bucket_id, , drop = FALSE]
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
        (norm(data, type = "F")^2 - colSums(data)^2 / n) / 2 / data_all_var +
          n / 2 * (log(data_all_var) + log(2 * pi))
      }
      mean_loss_result <- fastcpd(
        formula = ~ . - 1,
        data = data.frame(data),
        beta = (kDimension + 1) * log(nrow(data)) / 2,
        p = kDimension,
        cost = mean_loss
      )

      for (cp in mean_loss_result@cp_set) {
        ordinal_mapped_cp <-
          which(cumsum(cookie_bucket_id_list != cookie_bucket_id) == cp)[1]
        change_locations_cookie_bucket[[experiment_id]][[cookie_bucket_id]] <-
          ordinal_mapped_cp
      }
    }

    cp_for_eid <- Reduce("c", change_locations_cookie_bucket[[experiment_id]])
    d_capital <- mean(cp_for_eid)
    d_j <- rep(list(NULL), 20)
    for (j in seq_len(20)) {
      for (cookie_bucket_id in seq_len(20)) {
        if (j != cookie_bucket_id) {
          d_j[[j]] <- c(
            d_j[[j]],
            change_locations_cookie_bucket[[experiment_id]][[cookie_bucket_id]]
          )
        }
      }
    }
    d_j_bar <- sapply(d_j, mean)
    d_capital_j <- 20 * d_capital - 19 * d_j_bar
    d_bar <- mean(d_capital_j)
    sd_2 <- sum((d_capital_j - d_bar)^2) / 19
    ci <- c(d_capital - 2.093 * sqrt(sd_2) / sqrt(20), d_capital +
              2.093 * sqrt(sd_2) / sqrt(20))

    if (
      floor(ci[1]) <= 121 && ceiling(ci[2]) >= 121
    ) {
      containing_change_point[experiment_id] <- TRUE
    }
  }

  testthat::expect_equal(sum(containing_change_point), 927)
})

testthat::test_that("all examples in the documentation", {
  testthat::skip("This test is intended to be run manually.")
  fastcpd_documentation <- readLines("R/fastcpd.R")
  examples_index_start <-
    which(fastcpd_documentation == "#' # Linear regression")
  examples_index_end <-
    which(fastcpd_documentation == "#' summary(huber_regression_result)")
  examples <- fastcpd_documentation[examples_index_start:examples_index_end]
  for (i in seq_along(examples)) {
    if (fastcpd_documentation[i] == "#'") {
      fastcpd_documentation[i] <- ""
    } else {
      fastcpd_documentation[i] <-
        substr(fastcpd_documentation[i], 4, nchar(fastcpd_documentation[i]))
    }
  }
  source(textConnection(paste(
    fastcpd_documentation[examples_index_start:examples_index_end],
    collapse = "\n"
  )))
})
