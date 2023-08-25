test_that("linear regression", {
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

  expect_equal(result@cp_set, c(98, 202))
})

test_that("logistic regression", {
  # This is the same example with `fastcpd` documentation. Please keep it in
  # sync if the documentation ever changes.
  set.seed(1)

  kChangePointLocation <- 125
  kNumberOfDataPoints <- 300
  kDimension <- 5

  # There are 300 five-dimensional data points.
  x <- matrix(rnorm(kNumberOfDataPoints * kDimension, 0, 1), ncol = kDimension)

  # Randomly generate coefficients with different means.
  theta <- rbind(rnorm(kDimension, 0, 1), rnorm(kDimension, 2, 1))

  # Randomly generate response variables based on the segmented data and
  # corresponding coefficients
  y <- c(
    rbinom(kChangePointLocation, 1, 1 / (1 + exp(-x[1:kChangePointLocation, ] %*% theta[1, ]))),
    rbinom(kNumberOfDataPoints - kChangePointLocation, 1, 1 / (1 + exp(-x[(kChangePointLocation + 1):kNumberOfDataPoints, ] %*% theta[2, ])))
  )

  change_points_binomial_fastcpd <- suppressWarnings(fastcpd(
    formula = y ~ . - 1,
    data = data.frame(y = y, x = x),
    family = "binomial",
    segment_count = 5
  ))@cp_set

  expect_equal(change_points_binomial_fastcpd, kChangePointLocation)
})

test_that("poisson regression", {
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

  expect_equal(result@cp_set, c(329, 728, 1021, 1107, 1325))
})

test_that("penalized linear regression", {
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

  expect_equal(result@cp_set, c(300, 700, 1000))
})

test_that("confidence interval experiment", {
  set.seed(1)
  kDimension <- 1
  change_point_locations <- NULL
  for (experiment_id in seq_len(20)) {
    data <- rbind(
      mvtnorm::rmvnorm(300, mean = rep(0, kDimension), sigma = diag(100, kDimension)),
      mvtnorm::rmvnorm(400, mean = rep(50, kDimension), sigma = diag(100, kDimension)),
      mvtnorm::rmvnorm(300, mean = rep(2, kDimension), sigma = diag(100, kDimension))
    )
    segment_count_guess <- 10
    block_size <- max(floor(sqrt(nrow(data)) / (segment_count_guess + 1)), 2)
    block_count <- floor(nrow(data) / block_size)
    data_all_vars <- rep(0, block_count)
    for (block_index in seq_len(block_count)) {
      block_start <- (block_index - 1) * block_size + 1
      block_end <- if (block_index < block_count) block_index * block_size else nrow(data)
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

  change_point_locations_cookie_bucket <- NULL
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
      block_end <- if (block_index < block_count) block_index * block_size else nrow(data)
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
      ordinal_mapped_cp <- which(cumsum(cookie_bucket_id_list != cookie_bucket_id) == cp)[1]
      change_point_locations_cookie_bucket <- c(change_point_locations_cookie_bucket, ordinal_mapped_cp)
    }
  }
  
  table_change_point_locations <- table(change_point_locations)
  expect_equal(rownames(table_change_point_locations), c("269", "299", "300", "700"))
  expect_equal(unname(table_change_point_locations), c(1, 1, 19, 20), ignore_attr = TRUE)

  table_change_point_locations_cookie_bucket <- table(change_point_locations_cookie_bucket)
  expect_equal(rownames(table_change_point_locations_cookie_bucket), c("299", "300", "697", "700"))
  expect_equal(unname(table_change_point_locations_cookie_bucket), c(1, 19, 1, 19), ignore_attr = TRUE)
})
