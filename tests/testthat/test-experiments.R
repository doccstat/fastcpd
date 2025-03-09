testthat::test_that(
  "1d mean custom", {
    set.seed(1)
    result_mean <- fastcpd(
      formula = ~ . - 1,
      data = data.frame(x = c(rnorm(100, 0, 10), rnorm(100, 100, 10))),
      cost = function(data) {
        n <- nrow(data)
        demeaned_data <- sweep(data, 2, colMeans(data))
        n / 2 * (log(100) + log(2 * pi) + norm(demeaned_data, "2")^2 / n / 100)
      },
      cp_only = TRUE
    )
    testthat::expect_equal(result_mean@cp_set, 100)
  }
)

testthat::test_that(
  "warm start", {
    set.seed(1)
    p <- 1
    x <- matrix(rnorm(300 * p, 0, 1), ncol = p)
    theta <- rbind(rnorm(p, 0, 1), rnorm(p, 4, 1))
    y <- c(
      rbinom(150, 1, 1 / (1 + exp(-x[1:150, ] * theta[1, ]))),
      rbinom(150, 1, 1 / (1 + exp(-x[151:300, ] * theta[2, ])))
    )
    result <- suppressWarnings(
      fastcpd.binomial(
        data.frame(y = y, x = x),
        beta = "BIC",
        vanilla_percentage = 1,
        warm_start = TRUE,
        cost_adjustment = NULL
      )
    )
    testthat::expect_equal(result@cp_set, 134)
  }
)

testthat::test_that(
  "variance estimation", {
    set.seed(1)
    testthat::expect_message(
      try(fastcpd.ar(c(0, 0, 0, 0, 0, 0), 2), silent = TRUE),
      "Variance estimation failed for block 1."
    )
  }
)

testthat::skip("These tests are intended to be run manually.")

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
      multiple_epochs = function(segment_length) 1,
      epsilon = 1e-5,
      cost = logistic_loss,
      cost_gradient = logistic_loss_gradient,
      cost_hessian = logistic_loss_hessian
    )

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

    testthat::expect_equal(result_builtin@cp_set, 200)
    testthat::expect_equal(result_custom@cp_set, 201)
    testthat::expect_equal(result_custom_two_epochs@cp_set, 200)
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
      beta = (1 + 1 + 1) * log(n) / 2,
      p = 1 + 1,
      cost = ar1_loss
    )

    testthat::expect_equal(result@cp_set, 614)
  }
)

testthat::test_that(
  "confidence interval experiment", {
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
      mean_loss_result <- fastcpd::fastcpd(
        formula = ~ . - 1,
        data = data.frame(data),
        beta = (kDimension + 1) * log(nrow(data)) / 2,
        p = kDimension,
        cost = mean_loss
      )
      change_point_locations <-
        c(change_point_locations, mean_loss_result@cp_set)
    }

    cps_cookie_bucket <- NULL
    cookie_bucket_id_list <- sample.int(n = 20, size = 1000, replace = TRUE)
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
      mean_loss_result <- fastcpd::fastcpd(
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
        cps_cookie_bucket <-
          c(cps_cookie_bucket, ordinal_mapped_cp)
      }
    }

    table_change_point_locations <- table(change_point_locations)
    testthat::expect_equal(
      rownames(table_change_point_locations), c("269", "299", "300", "700")
    )
    testthat::expect_equal(
      unname(table_change_point_locations), c(1, 1, 19, 20), ignore_attr = TRUE
    )

    table_cp_cookie_bucket <- table(cps_cookie_bucket)
    testthat::expect_equal(
      rownames(table_cp_cookie_bucket), c("299", "300", "697", "700")
    )
    testthat::expect_equal(
      unname(table_cp_cookie_bucket), c(1, 19, 1, 19), ignore_attr = TRUE
    )
  }
)

testthat::test_that(  # nolint: cyclomatic complexity
  "confidence interval experiment with one change point", {
    set.seed(1)
    kDimension <- 1  # nolint: Google Style Guide
    change_point_locations <- list()
    cps_cookie_bucket <- rep(list(rep(list(NULL), 20)), 500)
    containing_change_point <- matrix(NA, 500, 5)

    # 500 experiments and count the number of times change point will be inside
    # the confidence interval to verify the logics of confidence interval.
    for (experiment_id in seq_len(500)) {
      data <- rbind(
        mvtnorm::rmvnorm(
          130, mean = rep(0, kDimension), sigma = diag(100, kDimension)
        ),
        mvtnorm::rmvnorm(
          130, mean = rep(50, kDimension), sigma = diag(100, kDimension)
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
        n / 2 * (
          log(data_all_var) + log(2 * pi) +
            sum((data - colMeans(data))^2 / data_all_var) / n
        )
      }
      mean_loss_result <- fastcpd::fastcpd(
        formula = ~ . - 1,
        data = data.frame(data),
        beta = (kDimension + 1) * log(nrow(data)) / 2,
        p = kDimension,
        cost = mean_loss
      )

      # Store the change point locations for each experiment as a baseline for
      # cookie bucket experiments.
      change_point_locations[[experiment_id]] <- mean_loss_result@cp_set

      # Randomly assign 20 cookie buckets to each experiment.
      cookie_bucket_id_list <-
        sample.int(n = 20, size = 130 + 130, replace = TRUE)
      all_data <- data

      # Do the cookie bucket experiment for each cookie bucket.
      for (cookie_bucket_id in seq_len(20)) {
        # Exclude the data from the cookie bucket.
        data <-
          all_data[cookie_bucket_id_list != cookie_bucket_id, , drop = FALSE]
        segment_count_guess <- 10
        block_size <-
          max(floor(sqrt(nrow(data)) / (segment_count_guess + 1)), 2)
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

        # Obtain the change points for the cookie bucket experiment.
        mean_loss_result <- fastcpd::fastcpd(
          formula = ~ . - 1,
          data = data.frame(data),
          beta = (kDimension + 1) * log(nrow(data)) / 2,
          p = kDimension,
          cost = mean_loss
        )

        # Map the change points to the original index.
        for (cp in mean_loss_result@cp_set) {
          ordinal_mapped_cp <-
            which(cumsum(cookie_bucket_id_list != cookie_bucket_id) == cp)[1]
          cps_cookie_bucket[[experiment_id]][[cookie_bucket_id]] <-
            ordinal_mapped_cp
        }
      }

      # Concatenate the change points from all cookie bucket experiments.
      cp_for_eid <- Reduce("c", cps_cookie_bucket[[experiment_id]])

      # Calculate the mean as the center of confidence interval.
      d_capital <- mean(cp_for_eid)

      # Concatenate the change ponints from all cookie bucket experiments
      # except the current one.
      d_j <- rep(list(NULL), 20)
      for (j in seq_len(20)) {
        for (cookie_bucket_id in seq_len(20)) {
          if (j != cookie_bucket_id) {
            d_j[[j]] <- c(
              d_j[[j]],
              cps_cookie_bucket[[experiment_id]][[cookie_bucket_id]]
            )
          }
        }
      }
      d_j_bar <- sapply(d_j, mean)

      # Calculate the jackknife cookie bucket change point.
      d_capital_j <- 20 * d_capital - 19 * d_j_bar
      d_bar <- mean(d_capital_j)
      sd_2 <- sum((d_capital_j - d_bar)^2) / 19

      # Following t-distribution with 19 degrees of freedom.
      containing_change_point[experiment_id, ] <- c(
        # 70% confidence interval.
        ceiling(d_capital - 1.066 * sqrt(sd_2) / sqrt(20)) <= 130 &&
          floor(d_capital + 1.066 * sqrt(sd_2) / sqrt(20)) >= 130,
        # 80% confidence interval.
        ceiling(d_capital - 1.328 * sqrt(sd_2) / sqrt(20)) <= 130 &&
          floor(d_capital + 1.328 * sqrt(sd_2) / sqrt(20)) >= 130,
        # 90% confidence interval.
        ceiling(d_capital - 1.729 * sqrt(sd_2) / sqrt(20)) <= 130 &&
          floor(d_capital + 1.729 * sqrt(sd_2) / sqrt(20)) >= 130,
        # 95% confidence interval.
        ceiling(d_capital - 2.093 * sqrt(sd_2) / sqrt(20)) <= 130 &&
          floor(d_capital + 2.093 * sqrt(sd_2) / sqrt(20)) >= 130,
        # 98% confidence interval.
        ceiling(d_capital - 2.539 * sqrt(sd_2) / sqrt(20)) <= 130 &&
          floor(d_capital + 2.539 * sqrt(sd_2) / sqrt(20)) >= 130
      )
    }

    testthat::expect_equal(
      colSums(containing_change_point), c(433, 436, 442, 449, 460)
    )
  }
)

testthat::test_that(  # nolint: cyclomatic complexity
  "confidence interval experiment with one change point for linear model", {
    set.seed(1)
    kDimension <- 3  # nolint: Google Style Guide
    change_point_locations <- list()
    cps_cookie_bucket <- rep(list(rep(list(NULL), 20)), 500)
    containing_change_point <- matrix(NA, 500, 5)

    # 500 experiments and count the number of times change point will be inside
    # the confidence interval to verify the logics of confidence interval.
    for (experiment_id in seq_len(500)) {
      x <- mvtnorm::rmvnorm(800, rep(0, p), diag(p))
      theta_0 <- rbind(c(1, 1.2, -1), c(-1, 0, 0.5))
      y <- c(
        x[1:500, ] %*% theta_0[1, ] + rnorm(100, 0, 1),
        x[501:800, ] %*% theta_0[2, ] + rnorm(100, 0, 1)
      )
      data <- data.frame(y = y, x = x)
      result <- fastcpd::fastcpd(
        formula = y ~ . - 1,
        data = data,
        family = "lm"
      )

      # Store the change point locations for each experiment as a baseline for
      # cookie bucket experiments.
      change_point_locations[[experiment_id]] <- result@cp_set

      # Randomly assign 20 cookie buckets to each experiment.
      cookie_bucket_id_list <-
        sample.int(n = 20, size = 800, replace = TRUE)
      all_data <- data

      # Do the cookie bucket experiment for each cookie bucket.
      for (cookie_bucket_id in seq_len(20)) {
        # Exclude the data from the cookie bucket.
        data <-
          all_data[cookie_bucket_id_list != cookie_bucket_id, , drop = FALSE]

        # Obtain the change points for the cookie bucket experiment.
        result <- fastcpd::fastcpd(
          formula = y ~ . - 1,
          data = data,
          family = "lm"
        )

        # Map the change points to the original index.
        for (cp in result@cp_set) {
          ordinal_mapped_cp <-
            which(cumsum(cookie_bucket_id_list != cookie_bucket_id) == cp)[1]
          cps_cookie_bucket[[experiment_id]][[cookie_bucket_id]] <-
            ordinal_mapped_cp
        }
      }

      # Concatenate the change points from all cookie bucket experiments.
      cp_for_eid <- Reduce("c", cps_cookie_bucket[[experiment_id]])

      # Calculate the mean as the center of confidence interval.
      d_capital <- mean(cp_for_eid)

      # Concatenate the change ponints from all cookie bucket experiments
      # except the current one.
      d_j <- rep(list(NULL), 20)
      for (j in seq_len(20)) {
        for (cookie_bucket_id in seq_len(20)) {
          if (j != cookie_bucket_id) {
            d_j[[j]] <- c(
              d_j[[j]],
              cps_cookie_bucket[[experiment_id]][[cookie_bucket_id]]
            )
          }
        }
      }
      d_j_bar <- sapply(d_j, mean)

      # Calculate the jackknife cookie bucket change point.
      d_capital_j <- 20 * d_capital - 19 * d_j_bar
      d_bar <- mean(d_capital_j)
      sd_2 <- sum((d_capital_j - d_bar)^2) / 19

      # Following t-distribution with 19 degrees of freedom.
      containing_change_point[experiment_id, ] <- c(
        # 70% confidence interval.
        ceiling(d_capital + 1.066 * sqrt(sd_2) / sqrt(20)) >= 500 &&
          floor(d_capital - 1.066 * sqrt(sd_2) / sqrt(20)) <= 500,
        # 80% confidence interval.
        ceiling(d_capital + 1.328 * sqrt(sd_2) / sqrt(20)) >= 500 &&
          floor(d_capital - 1.328 * sqrt(sd_2) / sqrt(20)) <= 500,
        # 90% confidence interval.
        ceiling(d_capital + 1.729 * sqrt(sd_2) / sqrt(20)) >= 500 &&
          floor(d_capital - 1.729 * sqrt(sd_2) / sqrt(20)) <= 500,
        # 95% confidence interval.
        ceiling(d_capital + 2.093 * sqrt(sd_2) / sqrt(20)) >= 500 &&
          floor(d_capital - 2.093 * sqrt(sd_2) / sqrt(20)) <= 500,
        # 98% confidence interval.
        ceiling(d_capital + 2.539 * sqrt(sd_2) / sqrt(20)) >= 500 &&
          floor(d_capital - 2.539 * sqrt(sd_2) / sqrt(20)) <= 500
      )
    }

    testthat::expect_equal(
      colSums(containing_change_point), c(371, 379, 390, 393, 399)
    )
  }
)
