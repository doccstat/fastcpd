testthat::test_that(
  "confidence interval experiment", {
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

testthat::test_that(
  "confidence interval experiment with one change point", {
    testthat::skip("This test is intended to be run manually.")
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

testthat::test_that(
  "confidence interval experiment with one change point for linear model", {
    testthat::skip("This test is intended to be run manually.")
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
        family = "gaussian"
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
          family = "gaussian"
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

testthat::test_that(
  "all examples in the documentation", {
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
  }
)

testthat::test_that(
  "build-in binomial performance on large data set with n = 10^4, p = 5", {
    testthat::skip("This test is intended to be run manually.")
    set.seed(1)
    n <- 10^4
    p <- 5
    segment_count <- n / 500
    theta_mean <- 5
    x <- mvtnorm::rmvnorm(n, mean = rep(0, p), sigma = diag(1, p))
    theta <- matrix(NA, segment_count, p)
    for (segment_count_index in seq_len(segment_count)) {
      theta[segment_count_index, ] <- rnorm(p, theta_mean, 5)
      theta_mean <- -theta_mean
    }
    y <- matrix(NA, n, 1)
    for (segment_count_index in seq_len(segment_count)) {
      segment_index <- (segment_count_index - 1) * 500 + seq_len(500)
      segment_theta <- theta[segment_count_index, ]
      y[segment_index] <-
        rbinom(500, 1, 1 / (1 + exp(-x[segment_index, ] %*% segment_theta)))
    }

    warning_messages <- testthat::capture_warnings(
      runtime <- system.time(
        result <- fastcpd::fastcpd(
          formula = y ~ . - 1,
          data = data.frame(y = y, x = x),
          family = "binomial"
        )
      )
    )

    testthat::expect_equal(
      warning_messages,
      rep("fit_glm: fitted probabilities numerically 0 or 1 occurred", 16)
    )

    # Discard the runtime value since it is different depending on the machine.
    # The run time is less than 10 seconds on a GitHub codespace with 2 cores.
    invisible(runtime)

    true_change_points <- 500 * seq_len(segment_count - 1)
    testthat::expect_equal(
      result@cp_set,
      true_change_points +
        c(1, 1, 0, 5, 0, 2, 0, -1, 2, 2, 1, 1, 2, 2, -4, 2, 10, 2, 0)
    )
  }
)

testthat::test_that(
  "build-in binomial performance on large data set with n = 10^4, p = 10", {
    testthat::skip("This test is intended to be run manually.")
    set.seed(1)
    n <- 10^4
    p <- 10
    segment_count <- n / 500
    theta_mean <- 5
    x <- mvtnorm::rmvnorm(n, mean = rep(0, p), sigma = diag(1, p))
    theta <- matrix(NA, segment_count, p)
    for (segment_count_index in seq_len(segment_count)) {
      theta[segment_count_index, ] <- rnorm(p, theta_mean, 5)
      theta_mean <- -theta_mean
    }
    y <- matrix(NA, n, 1)
    for (segment_count_index in seq_len(segment_count)) {
      segment_index <- (segment_count_index - 1) * 500 + seq_len(500)
      segment_theta <- theta[segment_count_index, ]
      y[segment_index] <-
        rbinom(500, 1, 1 / (1 + exp(-x[segment_index, ] %*% segment_theta)))
    }

    warning_messages <- testthat::capture_warnings(
      runtime <- system.time(
        result <- fastcpd::fastcpd(
          formula = y ~ . - 1,
          data = data.frame(y = y, x = x),
          family = "binomial"
        )
      )
    )

    testthat::expect_equal(
      warning_messages,
      rep("fit_glm: fitted probabilities numerically 0 or 1 occurred", 14)
    )

    # Discard the runtime value since it is different depending on the machine.
    # The run time is less than 10 seconds on a GitHub codespace with 2 cores.
    invisible(runtime)
    #    user  system elapsed
    #  10.444   0.106   7.765

    true_change_points <- 500 * seq_len(segment_count - 1)
    testthat::expect_equal(
      result@cp_set,
      true_change_points +
        c(5, 6, 5, 6, 5, 4, -4, 4, -1, 2, 1, 0, 2, 1, 5, 1, 4, 9, -1)
    )
  }
)

testthat::test_that(
  "build-in binomial performance on large data set with n = 10^4, p = 20", {
    testthat::skip("This test is intended to be run manually.")
    set.seed(1)
    n <- 10^4
    p <- 20
    segment_count <- n / 500
    theta_mean <- 5
    x <- mvtnorm::rmvnorm(n, mean = rep(0, p), sigma = diag(1, p))
    theta <- matrix(NA, segment_count, p)
    for (segment_count_index in seq_len(segment_count)) {
      theta[segment_count_index, ] <- rnorm(p, theta_mean, 5)
      theta_mean <- -theta_mean
    }
    y <- matrix(NA, n, 1)
    for (segment_count_index in seq_len(segment_count)) {
      segment_index <- (segment_count_index - 1) * 500 + seq_len(500)
      segment_theta <- theta[segment_count_index, ]
      y[segment_index] <-
        rbinom(500, 1, 1 / (1 + exp(-x[segment_index, ] %*% segment_theta)))
    }

    warning_messages <- testthat::capture_warnings(
      runtime <- system.time(
        result <- fastcpd::fastcpd(
          formula = y ~ . - 1,
          data = data.frame(y = y, x = x),
          family = "binomial"
        )
      )
    )

    testthat::expect_equal(
      warning_messages,
      rep("fit_glm: fitted probabilities numerically 0 or 1 occurred", 7)
    )

    # Discard the runtime value since it is different depending on the machine.
    # The run time is less than 10 seconds on a GitHub codespace with 2 cores.
    invisible(runtime)
    #   user  system elapsed
    # 36.672  14.923  28.899

    true_change_points <- 500 * seq_len(segment_count - 1)
    testthat::expect_equal(
      result@cp_set,
      true_change_points +
        c(11, -1, 8, 9, 12, 5, 9, -16, 28, 0, 0, 66, 26, 38, 7, 36, 8, 9, -3)
    )
  }
)

testthat::test_that(
  "build-in binomial performance on large data set with n = 3 * 10^4, p = 30", {
    testthat::skip("This test is intended to be run manually.")
    set.seed(1)
    n <- 3 * 10^4
    p <- 30
    segment_count <- n / 1000
    theta_mean <- 5
    x <- mvtnorm::rmvnorm(n, mean = rep(0, p), sigma = diag(1, p))
    theta <- matrix(NA, segment_count, p)
    for (segment_count_index in seq_len(segment_count)) {
      theta[segment_count_index, ] <- rnorm(p, theta_mean, 5)
      theta_mean <- -theta_mean
    }
    y <- matrix(NA, n, 1)
    for (segment_count_index in seq_len(segment_count)) {
      segment_index <- (segment_count_index - 1) * 1000 + seq_len(1000)
      segment_theta <- theta[segment_count_index, ]
      y[segment_index] <-
        rbinom(1000, 1, 1 / (1 + exp(-x[segment_index, ] %*% segment_theta)))
    }

    warning_messages <- testthat::capture_warnings(
      runtime <- system.time(
        result <- fastcpd::fastcpd(
          formula = y ~ . - 1,
          data = data.frame(y = y, x = x),
          family = "binomial"
        )
      )
    )

    testthat::expect_equal(
      warning_messages,
      rep("fit_glm: fitted probabilities numerically 0 or 1 occurred", 8)
    )

    # Discard the runtime value since it is different depending on the machine.
    # The run time is less than 10 seconds on a GitHub codespace with 2 cores.
    invisible(runtime)
    #    user  system elapsed
    # 683.446 356.068 576.744

    true_change_points <- 1000 * seq_len(segment_count - 1)
    testthat::expect_equal(
      result@cp_set,
      na.exclude(true_change_points + c(
        21, 50, 3, 91, 12, 27, NA, -317, 36, 168, 19, 79, 43, 97,
        0, 18, 5, 9, 63, 16, 26, 95, 4, -4, 24, 30, 57, 19, 19
      )),
      ignore_attr = TRUE
    )
  }
)

# testthat::test_that(
#   "real data analysis FRED macroeconomic", {
#     testthat::skip("This test is intended to be run manually.")
#     covariate_names <- c(
#       "UNRATE",
#       "PCEPILFE",
#       "HOUST",
#       "PAYEMS",
#       "RSXFS",
#       "INDPRO",
#       "FEDFUNDS"
#     )

#     # unemployment rate: UNRATE
#     # core PCE: PCEPILFE
#     # Housing starts: HOUST
#     # nonfarm employment: PAYEMS
#     # retail sales: RSXFS
#     # industrial production index: INDPRO
#     # Federal Funds rate: FEDFUNDS

#     covariates_url <- paste(
#       "https://fred.stlouisfed.org/graph/fredgraph.csv?id=",
#       covariate_names,
#       sep = ""
#     )

#     covariate_tables <- list()
#     min_date <- as.Date("1900-01-01")
#     max_date <- as.Date("2023-09-23")

#     for (covariate_name_index in seq_along(covariate_names)) {
#       covariate_tables[[covariate_names[covariate_name_index]]] <- read.csv(
#         covariates_url[covariate_name_index],
#         header = TRUE
#       )
#       min_date <- max(
#         min_date,
#         min(as.Date(
#           covariate_tables[[covariate_names[covariate_name_index]]]$DATE
#         ))
#       )
#       max_date <- min(
#         max_date,
#         max(as.Date(
#           covariate_tables[[covariate_names[covariate_name_index]]]$DATE
#         ))
#       )
#     }

#     covariate_table <- data.frame(date = seq(min_date, max_date, by = "months"))

#     for (covariate_name in covariate_names) {
#       covariate_table[[covariate_name]] <- covariate_tables[[covariate_name]][
#         covariate_tables[[covariate_name]]$DATE >= min_date &
#           covariate_tables[[covariate_name]]$DATE <= max_date,
#         covariate_name
#       ]
#     }

#     molten_covariate_table <- reshape2::melt(
#       covariate_table,
#       id.vars = "date",
#       variable.name = "covariate_name",
#       value.name = "covariate_value"
#     )

#     houst <- 100 * diff(log(covariate_table$HOUST))
#     n <- length(houst) - 2
#     y <- houst[2 + seq_len(n)]
#     x <- cbind(
#       houst[1 + seq_len(n)],
#       houst[seq_len(n)]
#     )
#     # n <- nrow(covariate_table) - 2
#     # y <- covariate_table[2 + seq_len(n), c("HOUST")]
#     # x <- cbind(
#     #   covariate_table[1 + seq_len(n), c("HOUST")],
#     #   covariate_table[seq_len(n), c("HOUST")]
#     # )
#     p <- ncol(x) / 2
#     multi_response_linear_loss <- function(data) {
#       x <- data[, (ncol(data) - p + 1):ncol(data)]
#       y <- data[, 1:(ncol(data) - p)]

#       if (nrow(data) <= p) {
#         x_t_x <- diag(p)
#       } else {
#         x_t_x <- crossprod(x)
#       }

#       norm(y - x %*% solve(x_t_x, t(x)) %*% y, type = "F")^2 / 2
#     }

#     data <- data.frame(y = y, x = x)
#     names(data) <- c(paste0("y.", seq_len(p)), paste0("x.", seq_len(2 * p)))
#     formula <- paste(
#       paste("cbind(", paste(paste0("y.", seq_len(p)), collapse = ", "), ")"),
#       "~",
#       paste(paste0("x.", seq_len(2 * p)), collapse = " + "),
#       "- 1"
#     )

#     result <- fastcpd::fastcpd(
#       formula = formula,
#       data = data,
#       beta = (2 * 2 * p + 1) * log(n) / 2,
#       cost = multi_response_linear_loss,
#       p = 2 * p,
#       trim = 0
#     )

#     testthat::expect_equal(
#       result@cp_set,
#       c(337, 338, 339)
#     )
#   }
# )

# testthat::test_that(
#   "variance estimation", {
#     testthat::skip("This test is intended to be run manually.")
#     set.seed(1)
#     p <- 3
#     n <- 300
#     cp <- c(100, 200)
#     x <- mvtnorm::rmvnorm(n, rep(0, p), diag(p))
#     theta_0 <- rbind(c(1, 1.2, -1), c(-1, 0, 0.5), c(0.5, -0.3, 0.2))
#     y <- c(
#       x[1:cp[1], ] %*% theta_0[1, ] + rnorm(cp[1], 0, sd = 3),
#       x[(cp[1] + 1):cp[2], ] %*% theta_0[2, ] + rnorm(cp[2] - cp[1], 0, sd = 3),
#       x[(cp[2] + 1):n, ] %*% theta_0[3, ] + rnorm(n - cp[2], 0, sd = 3)
#     )

#     pairwise_regression_cv <- function(
#       bandwidth_adjustment_factor, y, x, n, m
#     ) {
#       bandwidth <- bandwidth_adjustment_factor[1]
#       adjustment_factor <- bandwidth_adjustment_factor[2]
#       disjoint_count <- n / m
#       pairwise_regression_estimator <- rep(NA, disjoint_count)
#       for (v_i in seq_len(disjoint_count)) {
#         exclude_index <- seq_len(m) + (v_i - 1) * m
#         y_difference_matrix <- outer(y[-exclude_index], y[-exclude_index], "-")
#         x_difference_matrix <-
#           outer(c(x)[-exclude_index], c(x)[-exclude_index], "-")
#         s_matrix <- y_difference_matrix^2 / 2
#         d_matrix <- x_difference_matrix^2

#         z_array <-
#           y_difference_matrix[upper.tri(y_difference_matrix, diag = FALSE)]
#         lower_z <- quantile(z_array, 0.25)
#         upper_z <- quantile(z_array, 0.75)
#         z_lb <- lower_z - adjustment_factor * IQR(z_array)
#         z_ub <- upper_z + adjustment_factor * IQR(z_array)

#         a_set <- which(
#           d_matrix <= bandwidth & y_difference_matrix >= z_lb &
#             y_difference_matrix <= z_ub,
#           arr.ind = TRUE
#         )
#         a_set_n <- nrow(a_set)

#         s_1 <- sum(d_matrix[a_set])
#         s_2 <- sum(d_matrix[a_set]^2)
#         pairwise_regression_estimator[v_i] <- sum(
#           ((s_2 - s_1 * d_matrix) * s_matrix)[a_set]
#         ) / (a_set_n * s_2 - s_1^2)
#       }

#       y_difference_matrix <- outer(y, y, "-")
#       x_difference_matrix <- outer(c(x), c(x), "-")
#       s_matrix <- y_difference_matrix^2 / 2
#       d_matrix <- x_difference_matrix^2

#       z_array <-
#         y_difference_matrix[upper.tri(y_difference_matrix, diag = FALSE)]
#       lower_z <- quantile(z_array, 0.25)
#       upper_z <- quantile(z_array, 0.75)
#       z_lb <- lower_z - adjustment_factor * IQR(z_array)
#       z_ub <- upper_z + adjustment_factor * IQR(z_array)

#       a_set <- which(
#         d_matrix <= bandwidth & y_difference_matrix >= z_lb &
#           y_difference_matrix <= z_ub,
#         arr.ind = TRUE
#       )
#       a_set_n <- nrow(a_set)

#       s_1 <- sum(d_matrix[a_set])
#       s_2 <- sum(d_matrix[a_set]^2)
#       estimator_cv <- sum(
#         ((s_2 - s_1 * d_matrix) * s_matrix)[a_set]
#       ) / (a_set_n * s_2 - s_1^2)
#       sum((estimator_cv - pairwise_regression_estimator)^2)
#     }

#     baf_grid <- expand.grid(
#       4:8,
#       seq(2.0, 4.0, 0.4)
#     )

#     clusters <- parallel::makeCluster(parallel::detectCores() - 1)
#     doParallel::registerDoParallel(clusters)
#     baf_grid_values <- foreach::`%dopar%`(
#       foreach::foreach(baf_index = seq_len(nrow(baf_grid)), .combine = "c"), {
#         optim(
#           baf_grid[baf_index, ],
#           pairwise_regression_cv,
#           y = y,
#           x = x,
#           n = n,
#           m = 5
#         )$value
#       }
#     )
#     parallel::stopCluster(clusters)

#     baf_estimator <- baf_grid[which.min(baf_grid_values), ]
#     bandwidth <- baf_estimator$Var1
#     adjustment_factor <- baf_estimator$Var2

#     bandwidth <- 5
#     adjustment_factor <- 3

#     y_difference_matrix <- outer(y, y, "-")
#     x_difference_matrix <- as.matrix(dist(x))
#     s_matrix <- y_difference_matrix^2 / 2
#     d_matrix <- x_difference_matrix^2

#     z_array <-
#       y_difference_matrix[upper.tri(y_difference_matrix, diag = FALSE)]
#     lower_z <- quantile(z_array, 0.25)
#     upper_z <- quantile(z_array, 0.75)
#     z_lb <- lower_z - adjustment_factor * IQR(z_array)
#     z_ub <- upper_z + adjustment_factor * IQR(z_array)

#     a_set <- which(
#       d_matrix <= bandwidth & y_difference_matrix >= z_lb &
#         y_difference_matrix <= z_ub,
#       arr.ind = TRUE
#     )
#     a_set_n <- nrow(a_set)

#     s_1 <- sum(d_matrix[a_set])
#     s_2 <- sum(d_matrix[a_set]^2)
#     estimator_cv <-
#       sum(((s_2 - s_1 * d_matrix) * s_matrix)[a_set]) / (a_set_n * s_2 - s_1^2)

#     m <- 5
#     variance_estimation <- rep(NA, n - m)
#     for (i in 1:(n - m)) {
#       block_index <- seq_len(m) + i - 1
#       block_index_lagged <- seq_len(m) + i

#       y_block <- y[block_index]
#       x_block <- x[block_index, ]

#       y_block_lagged <- y[block_index_lagged]
#       x_block_lagged <- x[block_index_lagged, ]

#       x_t_x <- crossprod(x_block)
#       x_t_x_lagged <- crossprod(x_block_lagged)

#       block_slope <- solve(crossprod(x_block), crossprod(x_block, y_block))
#       block_lagged_slope <- solve(
#         crossprod(x_block_lagged), crossprod(x_block_lagged, y_block_lagged)
#       )
#       x_t_x_inv <- solve(x_t_x)
#       x_t_x_inv_lagged <- solve(x_t_x_lagged)
#       inv_product <- x_t_x_inv %*% x_t_x_inv_lagged
#       cross_term <-
#         inv_product %*% crossprod(x_block[-1, ], x_block_lagged[-m, ])
#       delta_numerator <- crossprod(block_slope - block_lagged_slope)
#       delta_denominator <-
#         sum(diag(x_t_x_inv + x_t_x_inv_lagged - 2 * cross_term))
#       variance_estimation[i] <- delta_numerator / delta_denominator
#     }

#     c(
#       sqrt(mean(variance_estimation)),
#       sqrt(estimator_cv)
#     )
#   }
# )

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
      formula = cbind(y.1, y.2) ~ x.1 + x.2 + x.3 + x.4 - 1,
      data = data.frame(y = data[, 1:2], x = data[, 3:6]),
      beta = (2 * 4 + 1) * log(n) / 2,
      cost = multi_response_linear_loss
    )

    testthat::expect_equal(result_custom@cp_set, 609)

    testthat::skip("TODO(doccstat): Fix VAR(p) model.")

    result_ts <- fastcpd_ts(x, "var", 2)

    testthat::expect_equal(result_ts@cp_set, c(424, 609, 720))

    result_var <- fastcpd(
      formula = ~ . - 1,
      data = data.frame(x = x),
      family = "var",
      order = 2
    )

    testthat::expect_equal(result_var@cp_set, 609)
  }
)
