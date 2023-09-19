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

# Everything in this script is provided as is. The purpose of this script is to
# do a sanity check on the C++ implementation of `fastcpd`.

# nolint start: script provided as is

#' Cost function for Logistic regression, i.e. binomial family in GLM.
#'
#' @param data Data to be used to calculate the cost values. The last column is
#'     the response variable.
#' @param family Family of the distribution.
#' @keywords internal
#'
#' @noRd
#' @return Cost value for the corresponding segment of data.
cost_glm_binomial <- function(data, family = "binomial") {
  data <- as.matrix(data)
  p <- dim(data)[2] - 1
  out <- fastglm::fastglm(
    as.matrix(data[, 1:p]), data[, p + 1],
    family = family
  )
  return(out$deviance / 2)
}

#' Implementation of vanilla PELT for logistic regression type data.
#'
#' @param data Data to be used for change point detection.
#' @param beta Penalty coefficient for the number of change points.
#' @param cost Cost function to be used to calculate cost values.
#' @keywords internal
#'
#' @noRd
#' @return A list consisting of two: change point locations and negative log
#'     likelihood values for each segment.
pelt_vanilla_binomial <- function(data, beta, cost = cost_glm_binomial) {
  n <- dim(data)[1]
  p <- dim(data)[2] - 1
  Fobj <- c(-beta, 0)
  cp_set <- list(NULL, 0)
  set <- c(0, 1)
  for (t in 2:n)
  {
    m <- length(set)
    cval <- rep(NA, m)
    for (i in 1:m)
    {
      k <- set[i] + 1
      cval[i] <- 0
      if (t - k >= p - 1) cval[i] <- suppressWarnings(cost(data[k:t, ]))
    }
    obj <- cval + Fobj[set + 1] + beta
    min_val <- min(obj)
    ind <- which(obj == min_val)[1]
    cp_set_add <- c(cp_set[[set[ind] + 1]], set[ind])
    cp_set <- append(cp_set, list(cp_set_add))
    ind2 <- (cval + Fobj[set + 1]) <= min_val
    set <- c(set[ind2], t)
    Fobj <- c(Fobj, min_val)
  }
  cp <- cp_set[[n + 1]]
  nLL <- 0
  cp_loc <- unique(c(0, cp, n))
  for (i in 1:(length(cp_loc) - 1))
  {
    seg <- (cp_loc[i] + 1):cp_loc[i + 1]
    data_seg <- data[seg, ]
    out <- fastglm::fastglm(
      as.matrix(data_seg[, 1:p]), data_seg[, p + 1],
      family = "binomial"
    )
    nLL <- out$deviance / 2 + nLL
  }

  output <- list(cp, nLL)
  names(output) <- c("cp", "nLL")
  return(output)
}

#' Function to update the coefficients using gradient descent.
#'
#' @param data_new New data point used to update the coeffient.
#' @param coef Previous coeffient to be updated.
#' @param cum_coef Summation of all the past coefficients to be used in
#'     averaging.
#' @param cmatrix Hessian matrix in gradient descent.
#' @param epsilon Small adjustment to avoid singularity when doing inverse on
#'     the Hessian matrix.
#' @keywords internal
#'
#' @noRd
#' @return A list of values containing the new coefficients, summation of
#'     coefficients so far and all the Hessian matrices.
cost_logistic_update <- function(
    data_new, coef, cum_coef, cmatrix, epsilon = 1e-10) {
  p <- length(data_new) - 1
  X_new <- data_new[1:p]
  Y_new <- data_new[p + 1]
  eta <- X_new %*% coef
  mu <- 1 / (1 + exp(-eta))
  cmatrix <- cmatrix + (X_new %o% X_new) * as.numeric((1 - mu) * mu)
  lik_dev <- as.numeric(-(Y_new - mu)) * X_new
  coef <- coef - solve(cmatrix + epsilon * diag(1, p), lik_dev)
  cum_coef <- cum_coef + coef
  return(list(coef, cum_coef, cmatrix))
}

#' Calculate negative log likelihood given data segment and guess of
#' coefficient.
#'
#' @param data Data to be used to calculate the negative log likelihood.
#' @param b Guess of the coefficient.
#' @keywords internal
#'
#' @noRd
#' @return Negative log likelihood.
neg_log_lik_binomial <- function(data, b) {
  p <- dim(data)[2] - 1
  X <- data[, 1:p, drop = FALSE]
  Y <- data[, p + 1, drop = FALSE]
  u <- as.numeric(X %*% b)
  L <- -Y * u + log(1 + exp(u))
  return(sum(L))
}

#' Find change points using dynamic programming with pruning and SeGD.
#'
#' @param data Data used to find change points.
#' @param beta Penalty coefficient for the number of change points.
#' @param B Initial guess on the number of change points.
#' @param trim Propotion of the data to ignore the change points at the
#'     beginning, ending and between change points.
#' @keywords internal
#'
#' @noRd
#' @return A list containing potential change point locations and negative log
#'     likelihood for each segment based on the change points guess.
segd_binomial <- function(data, beta, B = 10, trim = 0.025) {
  n <- dim(data)[1]
  p <- dim(data)[2] - 1
  Fobj <- c(-beta, 0)
  cp_set <- list(NULL, 0)
  set <- c(0, 1)

  # choose the initial values based on pre-segmentation

  index <- rep(1:B, rep(n / B, B))
  coef.int <- matrix(NA, B, p)
  for (i in 1:B)
  {
    out <- fastglm::fastglm(
      as.matrix(data[index == i, 1:p]),
      data[index == i, p + 1],
      family = "binomial"
    )
    coef.int[i, ] <- coef(out)
  }
  X1 <- data[1, 1:p]
  cum_coef <- coef <- matrix(coef.int[1, ], p, 1)
  e_eta <- exp(coef %*% X1)
  const <- e_eta / (1 + e_eta)^2
  cmatrix <- array((X1 %o% X1) * as.numeric(const), c(p, p, 1))

  for (t in 2:n)
  {
    m <- length(set)
    cval <- rep(NA, m)

    for (i in 1:(m - 1))
    {
      coef_c <- coef[, i]
      cum_coef_c <- cum_coef[, i]
      cmatrix_c <- cmatrix[, , i]
      out <- cost_logistic_update(data[t, ], coef_c, cum_coef_c, cmatrix_c)
      coef[, i] <- out[[1]]
      cum_coef[, i] <- out[[2]]
      cmatrix[, , i] <- out[[3]]
      k <- set[i] + 1
      cval[i] <- 0
      if (t - k >= p - 1) {
        cval[i] <-
          neg_log_lik_binomial(data[k:t, ], cum_coef[, i] / (t - k + 1))
      }
    }

    # the choice of initial values requires further investigation

    cval[m] <- 0
    Xt <- data[t, 1:p]
    cum_coef_add <- coef_add <- coef.int[index[t], ]
    e_eta_t <- exp(coef_add %*% Xt)
    const <- e_eta_t / (1 + e_eta_t)^2
    cmatrix_add <- (Xt %o% Xt) * as.numeric(const)

    coef <- cbind(coef, coef_add)
    cum_coef <- cbind(cum_coef, cum_coef_add)
    cmatrix <- abind::abind(cmatrix, cmatrix_add, along = 3)

    # Adding a momentum term (TBD)

    obj <- cval + Fobj[set + 1] + beta
    min_val <- min(obj)
    ind <- which(obj == min_val)[1]
    cp_set_add <- c(cp_set[[set[ind] + 1]], set[ind])
    cp_set <- append(cp_set, list(cp_set_add))
    ind2 <- (cval + Fobj[set + 1]) <= min_val
    set <- c(set[ind2], t)
    coef <- coef[, ind2, drop = FALSE]
    cum_coef <- cum_coef[, ind2, drop = FALSE]
    cmatrix <- cmatrix[, , ind2, drop = FALSE]
    Fobj <- c(Fobj, min_val)
  }

  # Remove change-points close to the boundaries

  cp <- cp_set[[n + 1]]
  if (length(cp) > 0) {
    ind3 <- (seq_len(length(cp)))[(cp < trim * n) | (cp > (1 - trim) * n)]
    cp <- cp[-ind3]
  }

  nLL <- 0
  cp_loc <- unique(c(0, cp, n))
  for (i in 1:(length(cp_loc) - 1))
  {
    seg <- (cp_loc[i] + 1):cp_loc[i + 1]
    data_seg <- data[seg, ]
    out <- fastglm::fastglm(
      as.matrix(data_seg[, 1:p]), data_seg[, p + 1],
      family = "binomial"
    )
    nLL <- out$deviance / 2 + nLL
  }

  output <- list(cp, nLL)
  names(output) <- c("cp", "nLL")
  return(output)
}

testthat::test_that("logistic regression", {
  testthat::skip("This test is intended to be run manually.")
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

  change_points_binomial_fastcpd <- suppressWarnings(fastcpd::fastcpd(
    formula = y ~ . - 1,
    data = data.frame(y = y, x = x),
    family = "binomial",
    segment_count = 5
  ))@cp_set

  change_points_binomial_fastcpd_sanity <- segd_binomial(
    cbind(x, y), (kDimension + 1) * log(kNumberOfDataPoints) / 2,
    B = 5
  )$cp

  testthat::expect_equal(
    change_points_binomial_fastcpd,
    kChangePointLocation
  )

  testthat::expect_equal(
    change_points_binomial_fastcpd,
    change_points_binomial_fastcpd_sanity
  )

  warning_messages <- testthat::capture_warnings(
    change_points_binomial_fastcpd_vanilla <- fastcpd::fastcpd(
      formula = y ~ . - 1,
      data = data.frame(y = y, x = x),
      family = "binomial",
      segment_count = 5,
      vanilla_percentage = 1
    )@cp_set
  )

  testthat::expect_equal(
    sort(warning_messages),
    rep(c(
      "fit_glm: algorithm did not converge",
      "fit_glm: fitted probabilities numerically 0 or 1 occurred"
    ), c(6, 3984))
  )

  change_points_binomial_fastcpd_vanilla_sanity <- pelt_vanilla_binomial(
    cbind(x, y), (kDimension + 1) * log(kNumberOfDataPoints) / 2
  )$cp

  testthat::expect_equal(
    change_points_binomial_fastcpd_vanilla,
    kChangePointLocation
  )
  testthat::expect_equal(
    c(0, change_points_binomial_fastcpd_vanilla),
    change_points_binomial_fastcpd_vanilla_sanity
  )
})

#' Cost function for Poisson regression.
#'
#' @param data Data to be used to calculate the cost values. The last column is
#'     the response variable.
#' @param family Family of the distribution.
#' @keywords internal
#'
#' @noRd
#' @return Cost value for the corresponding segment of data.
cost_glm_poisson <- function(data, family = "poisson") {
  data <- as.matrix(data)
  p <- dim(data)[2] - 1
  out <- fastglm(as.matrix(data[, 1:p]), data[, p + 1], family = family)
  return(out$deviance / 2)
}

#' Implementation of vanilla PELT for poisson regression type data.
#'
#' @param data Data to be used for change point detection.
#' @param beta Penalty coefficient for the number of change points.
#' @param cost Cost function to be used to calculate cost values.
#' @keywords internal
#'
#' @noRd
#' @return A list consisting of two: change point locations and negative log
#'     likelihood values for each segment.
pelt_vanilla_poisson <- function(data, beta, cost = cost_glm_poisson) {
  n <- dim(data)[1]
  p <- dim(data)[2] - 1
  Fobj <- c(-beta, 0)
  cp_set <- list(NULL, 0)
  set <- c(0, 1)
  for (t in 2:n)
  {
    m <- length(set)
    cval <- rep(NA, m)
    for (i in 1:m)
    {
      k <- set[i] + 1
      if (t - k >= p - 1) cval[i] <- suppressWarnings(cost(data[k:t, ])) else cval[i] <- 0
    }
    obj <- cval + Fobj[set + 1] + beta
    min_val <- min(obj)
    ind <- which(obj == min_val)[1]
    cp_set_add <- c(cp_set[[set[ind] + 1]], set[ind])
    cp_set <- append(cp_set, list(cp_set_add))
    ind2 <- (cval + Fobj[set + 1]) <= min_val
    set <- c(set[ind2], t)
    Fobj <- c(Fobj, min_val)
    # if (t %% 100 == 0) print(t)
  }
  cp <- cp_set[[n + 1]]
  output <- list(cp)
  names(output) <- c("cp")

  return(output)
}

#' Function to update the coefficients using gradient descent.
#'
#' @param data_new New data point used to update the coeffient.
#' @param coef Previous coeffient to be updated.
#' @param cum_coef Summation of all the past coefficients to be used in
#'     averaging.
#' @param cmatrix Hessian matrix in gradient descent.
#' @param epsilon Small adjustment to avoid singularity when doing inverse on
#'     the Hessian matrix.
#' @param G Upper bound for the coefficient.
#' @param L Winsorization lower bound.
#' @param H Winsorization upper bound.
#' @keywords internal
#'
#' @noRd
#' @return A list of values containing the new coefficients, summation of
#'     coefficients so far and all the Hessian matrices.
cost_poisson_update <- function(data_new, coef, cum_coef, cmatrix, epsilon = 0.001, G = 10^10, L = -20, H = 20) {
  p <- length(data_new) - 1
  X_new <- data_new[1:p]
  Y_new <- data_new[p + 1]
  eta <- X_new %*% coef
  mu <- exp(eta)
  cmatrix <- cmatrix + (X_new %o% X_new) * min(as.numeric(mu), G)
  lik_dev <- as.numeric(-(Y_new - mu)) * X_new
  coef <- coef - solve(cmatrix + epsilon * diag(1, p), lik_dev)
  coef <- Winsorize(coef, minval = L, maxval = H)
  cum_coef <- cum_coef + coef
  return(list(coef, cum_coef, cmatrix))
}

#' Calculate negative log likelihood given data segment and guess of
#' coefficient.
#'
#' @param data Data to be used to calculate the negative log likelihood.
#' @param b Guess of the coefficient.
#' @keywords internal
#'
#' @noRd
#' @return Negative log likelihood.
neg_log_lik_poisson <- function(data, b) {
  p <- dim(data)[2] - 1
  X <- data[, 1:p, drop = FALSE]
  Y <- data[, p + 1, drop = FALSE]
  u <- as.numeric(X %*% b)
  L <- -Y * u + exp(u) + lfactorial(Y)
  return(sum(L))
}

#' Find change points using dynamic programming with pruning and SeGD.
#'
#' @param data Data used to find change points.
#' @param beta Penalty coefficient for the number of change points.
#' @param B Initial guess on the number of change points.
#' @param trim Propotion of the data to ignore the change points at the
#'     beginning, ending and between change points.
#' @param epsilon Small adjustment to avoid singularity when doing inverse on
#'    the Hessian matrix.
#' @param G Upper bound for the coefficient.
#' @param L Winsorization lower bound.
#' @param H Winsorization upper bound.
#' @keywords internal
#'
#' @noRd
#' @return A list containing potential change point locations and negative log
#'     likelihood for each segment based on the change points guess.
segd_poisson <- function(data, beta, B = 10, trim = 0.03, epsilon = 0.001, G = 10^10, L = -20, H = 20) {
  n <- dim(data)[1]
  p <- dim(data)[2] - 1
  Fobj <- c(-beta, 0)
  cp_set <- list(NULL, 0)
  set <- c(0, 1)

  # choose the initial values based on pre-segmentation

  index <- rep(1:B, rep(n / B, B))
  coef.int <- matrix(NA, B, p)
  for (i in 1:B)
  {
    out <- fastglm(x = as.matrix(data[index == i, 1:p]), y = data[index == i, p + 1], family = "poisson")
    coef.int[i, ] <- coef(out)
  }
  X1 <- data[1, 1:p]
  cum_coef <- coef <- Winsorize(matrix(coef.int[1, ], p, 1), minval = L, maxval = H)
  e_eta <- exp(coef %*% X1)
  const <- e_eta
  cmatrix <- array((X1 %o% X1) * as.numeric(const), c(p, p, 1))

  for (t in 2:n)
  {
    m <- length(set)
    cval <- rep(NA, m)
    for (i in 1:(m - 1))
    {
      coef_c <- coef[, i]
      cum_coef_c <- cum_coef[, i]
      cmatrix_c <- cmatrix[, , i]
      out <- cost_poisson_update(data[t, ], coef_c, cum_coef_c, cmatrix_c, epsilon = epsilon, G = G, L = L, H = H)
      coef[, i] <- out[[1]]
      cum_coef[, i] <- out[[2]]
      cmatrix[, , i] <- out[[3]]
      k <- set[i] + 1
      cum_coef_win <- Winsorize(cum_coef[, i] / (t - k + 1), minval = L, maxval = H)
      if (t - k >= p - 1) cval[i] <- neg_log_lik_poisson(data[k:t, ], cum_coef_win) else cval[i] <- 0
    }

    # the choice of initial values requires further investigation

    cval[m] <- 0
    Xt <- data[t, 1:p]
    cum_coef_add <- coef_add <- Winsorize(coef.int[index[t], ], minval = L, maxval = H) ####
    e_eta_t <- exp(coef_add %*% Xt)
    const <- e_eta_t
    cmatrix_add <- (Xt %o% Xt) * as.numeric(const)
    coef <- cbind(coef, coef_add)
    cum_coef <- cbind(cum_coef, cum_coef_add)
    cmatrix <- abind::abind(cmatrix, cmatrix_add, along = 3)

    # Adding a momentum term (TBD)

    obj <- cval + Fobj[set + 1] + beta
    min_val <- min(obj)
    ind <- which(obj == min_val)[1]
    cp_set_add <- c(cp_set[[set[ind] + 1]], set[ind])
    cp_set <- append(cp_set, list(cp_set_add))
    ind2 <- (cval + Fobj[set + 1]) <= min_val
    set <- c(set[ind2], t)
    coef <- coef[, ind2, drop = FALSE]
    cum_coef <- cum_coef[, ind2, drop = FALSE]
    cmatrix <- cmatrix[, , ind2, drop = FALSE]
    Fobj <- c(Fobj, min_val)
  }

  # Remove change-points close to the boundaries

  cp <- cp_set[[n + 1]]
  if (length(cp) > 0) {
    ind3 <- (1:length(cp))[(cp < trim * n) | (cp > (1 - trim) * n)]
    if (length(ind3) > 0) cp <- cp[-ind3]
  }

  cp <- sort(unique(c(0, cp)))
  index <- which((diff(cp) < trim * n) == TRUE)
  if (length(index) > 0) cp <- floor((cp[-(index + 1)] + cp[-index]) / 2)
  cp <- cp[cp > 0]

  # nLL <- 0
  # cp_loc <- unique(c(0,cp,n))
  # for(i in 1:(length(cp_loc)-1))
  # {
  #   seg <- (cp_loc[i]+1):cp_loc[i+1]
  #   data_seg <- data[seg,]
  #   out <- fastglm(as.matrix(data_seg[, 1:p]), data_seg[, p+1], family="Poisson")
  #   nLL <- out$deviance/2 + nLL
  # }

  # output <- list(cp, nLL)
  # names(output) <- c("cp", "nLL")

  output <- list(cp)
  names(output) <- c("cp")

  return(output)
}

# Generate data from poisson regression models with change-points
#' @param n Number of observations.
#' @param d Dimension of the covariates.
#' @param true.coef True regression coefficients.
#' @param true.cp.loc True change-point locations.
#' @param Sigma Covariance matrix of the covariates.
#' @keywords internal
#'
#' @noRd
#' @return A list containing the generated data and the true cluster
#'    assignments.
data_gen_poisson <- function(n, d, true.coef, true.cp.loc, Sigma) {
  loc <- unique(c(0, true.cp.loc, n))
  if (dim(true.coef)[2] != length(loc) - 1) stop("true.coef and true.cp.loc do not match")
  x <- mvtnorm::rmvnorm(n, mean = rep(0, d), sigma = Sigma)
  y <- NULL
  for (i in 1:(length(loc) - 1))
  {
    mu <- exp(x[(loc[i] + 1):loc[i + 1], , drop = FALSE] %*% true.coef[, i, drop = FALSE])
    group <- rpois(length(mu), mu)
    y <- c(y, group)
  }
  data <- cbind(x, y)
  true_cluster <- rep(1:(length(loc) - 1), diff(loc))
  result <- list(data, true_cluster)
  return(result)
}

testthat::test_that("poisson regression", {
  testthat::skip("This test is intended to be run manually.")
  set.seed(1)
  n <- 1500
  d <- 5
  rho <- 0.9
  Sigma <- array(0, c(d, d))
  for (i in 1:d) {
    Sigma[i, ] <- rho^(abs(i - (1:d)))
  }
  delta <- c(5, 7, 9, 11, 13)
  a.sq <- 1
  delta.new <- delta * sqrt(a.sq) / sqrt(as.numeric(t(delta) %*% Sigma %*% delta))
  true.cp.loc <- c(375, 750, 1125)

  # regression coefficients

  true.coef <- matrix(0, nrow = d, ncol = length(true.cp.loc) + 1)
  true.coef[, 1] <- c(1, 1.2, -1, 0.5, -2)
  true.coef[, 2] <- true.coef[, 1] + delta.new
  true.coef[, 3] <- true.coef[, 1]
  true.coef[, 4] <- true.coef[, 3] - delta.new

  out <- data_gen_poisson(n, d, true.coef, true.cp.loc, Sigma)
  data <- out[[1]]
  g_tr <- out[[2]]
  beta <- log(n) * (d + 1) / 2

  change_points_poisson_fastcpd <- fastcpd::fastcpd(
    formula = y ~ . - 1,
    data = data.frame(y = data[, d + 1], x = data[, 1:d]),
    epsilon = 0.001,
    family = "poisson",
    segment_count = 10
  )@cp_set
  change_points_poisson_fastcpd_sanity <-
    segd_poisson(data, beta, trim = 0.03, B = 10, epsilon = 0.001, G = 10^10, L = -20, H = 20)$cp

  warning_messages <- testthat::capture_warnings(
    change_points_poisson_fastcpd_vanilla <- fastcpd::fastcpd(
      formula = y ~ . - 1,
      data = data.frame(y = data[, d + 1], x = data[, 1:d]),
      family = "poisson",
      segment_count = 10,
      vanilla_percentage = 1
    )@cp_set
  )

  testthat::expect_equal(
    warning_messages, rep("fit_glm: fitted rates numerically 0 occurred", 1539)
  )

  change_points_poisson_fastcpd_vanilla_sanity <-
    pelt_vanilla_poisson(data, beta)$cp

  testthat::expect_equal(
    change_points_poisson_fastcpd,
    c(380, 751, 1136, 1251)
  )

  testthat::expect_equal(
    change_points_poisson_fastcpd,
    change_points_poisson_fastcpd_sanity
  )

  testthat::expect_equal(
    change_points_poisson_fastcpd_vanilla,
    c(374, 752, 1133)
  )

  testthat::expect_equal(
    c(0, change_points_poisson_fastcpd_vanilla),
    change_points_poisson_fastcpd_vanilla_sanity
  )
})

#' Cost function for penalized linear regression.
#'
#' @param data Data to be used to calculate the cost values. The last column is
#'     the response variable.
#' @param lambda Penalty coefficient.
#' @param family Family of the distribution.
#' @keywords internal
#'
#' @noRd
#' @return Cost value for the corresponding segment of data.
cost_lasso <- function(data, lambda, family = "gaussian") {
  data <- as.matrix(data)
  n <- dim(data)[1]
  p <- dim(data)[2] - 1
  out <- glmnet(as.matrix(data[, 1:p]), data[, p + 1], family = family, lambda = lambda)
  return(deviance(out) / 2)
}

#' Implementation of vanilla PELT for penalized linear regression type data.
#'
#' @param data Data to be used for change point detection.
#' @param beta Penalty coefficient for the number of change points.
#' @param B Initial guess on the number of change points.
#' @param cost Cost function to be used to calculate cost values.
#' @param family Family of the distribution.
#' @keywords internal
#'
#' @noRd
#' @return A list consisting of two: change point locations and negative log
#'     likelihood values for each segment.
pelt_vanilla_lasso <- function(data, beta, B = 10, cost = cost_lasso, family = "gaussian") {
  n <- dim(data)[1]
  p <- dim(data)[2] - 1
  index <- rep(1:B, rep(n / B, B))
  err_sd <- act_num <- rep(NA, B)
  for (i in 1:B)
  {
    cvfit <- cv.glmnet(as.matrix(data[index == i, 1:p]), data[index == i, p + 1], family = family)
    coef <- coef(cvfit, s = "lambda.1se")[-1]
    resi <- data[index == i, p + 1] - as.matrix(data[index == i, 1:p]) %*% as.numeric(coef)
    err_sd[i] <- sqrt(mean(resi^2))
    act_num[i] <- sum(abs(coef) > 0)
  }
  err_sd_mean <- mean(err_sd) # only works if error sd is unchanged.
  act_num_mean <- mean(act_num)
  beta <- (act_num_mean + 1) * beta # seems to work but there might be better choices

  Fobj <- c(-beta, 0)
  cp_set <- list(NULL, 0)
  set <- c(0, 1)
  for (t in 2:n)
  {
    m <- length(set)
    cval <- rep(NA, m)
    for (i in 1:m)
    {
      k <- set[i] + 1
      if (t - k >= 1) cval[i] <- suppressWarnings(cost(data[k:t, ], lambda = err_sd_mean * sqrt(2 * log(p) / (t - k + 1)))) else cval[i] <- 0
    }
    obj <- cval + Fobj[set + 1] + beta
    min_val <- min(obj)
    ind <- which(obj == min_val)[1]
    cp_set_add <- c(cp_set[[set[ind] + 1]], set[ind])
    cp_set <- append(cp_set, list(cp_set_add))
    ind2 <- (cval + Fobj[set + 1]) <= min_val
    set <- c(set[ind2], t)
    Fobj <- c(Fobj, min_val)
    if (t %% 100 == 0) print(t)
  }
  cp <- cp_set[[n + 1]]
  # nLL <- 0
  # cp_loc <- unique(c(0,cp,n))
  # for(i in 1:(length(cp_loc)-1))
  # {
  #  seg <- (cp_loc[i]+1):cp_loc[i+1]
  #  data_seg <- data[seg,]
  #  out <- glmnet(as.matrix(data_seg[, 1:p]), data_seg[, p+1], lambda=lambda, family=family)
  #  nLL <- deviance(out)/2 + nLL
  # }
  # output <- list(cp, nLL)
  output <- list(cp)
  names(output) <- c("cp")
  return(output)
}

#' Function to update the coefficients using gradient descent.
#' @param a Coefficient to be updated.
#' @param lambda Penalty coefficient.
#' @keywords internal
#'
#' @noRd
#' @return Updated coefficient.
soft_threshold <- function(a, lambda) {
  sign(a) * pmax(abs(a) - lambda, 0)
}

#' Function to update the coefficients using gradient descent.
#'
#' @param data_new New data point used to update the coeffient.
#' @param coef Previous coeffient to be updated.
#' @param cum_coef Summation of all the past coefficients to be used in
#'     averaging.
#' @param cmatrix Hessian matrix in gradient descent.
#' @param lambda Penalty coefficient.
#' @keywords internal
#'
#' @noRd
#' @return A list of values containing the new coefficients, summation of
#'     coefficients so far and all the Hessian matrices.
cost_lasso_update <- function(data_new, coef, cum_coef, cmatrix, lambda) {
  p <- length(data_new) - 1
  X_new <- data_new[1:p]
  Y_new <- data_new[p + 1]
  mu <- X_new %*% coef
  cmatrix <- cmatrix + X_new %o% X_new
  # B <- as.vector(cmatrix_inv%*%X_new)
  # cmatrix_inv <- cmatrix_inv - B%o%B/(1+sum(X_new*B))
  lik_dev <- as.numeric(-(Y_new - mu)) * X_new
  coef <- coef - solve(cmatrix, lik_dev)
  nc <- norm(cmatrix, type = "F") # the choice of norm affects the speed. Spectral norm is more accurate but slower than F norm.
  coef <- soft_threshold(coef, lambda / nc)
  cum_coef <- cum_coef + coef
  return(list(coef, cum_coef, cmatrix))
}

#' Calculate negative log likelihood given data segment and guess of
#' coefficient.
#'
#' @param data Data to be used to calculate the negative log likelihood.
#' @param b Guess of the coefficient.
#' @param lambda Penalty coefficient.
#' @keywords internal
#'
#' @noRd
#' @return Negative log likelihood.
neg_log_lik_lasso <- function(data, b, lambda) {
  p <- dim(data)[2] - 1
  X <- data[, 1:p, drop = FALSE]
  Y <- data[, p + 1, drop = FALSE]
  resi <- Y - X %*% b
  L <- sum(resi^2) / 2 + lambda * sum(abs(b))
  return(L)
}

#' Find change points using dynamic programming with pruning and SeGD.
#'
#' @param data Data used to find change points.
#' @param beta Penalty coefficient for the number of change points.
#' @param B Initial guess on the number of change points.
#' @param trim Propotion of the data to ignore the change points at the
#'     beginning, ending and between change points.
#' @param epsilon Small adjustment to avoid singularity when doing inverse on
#'    the Hessian matrix.
#' @param family Family of the distribution.
#' @keywords internal
#'
#' @noRd
#' @return A list containing potential change point locations and negative log
#'     likelihood for each segment based on the change points guess.
segd_lasso <- function(data, beta, B = 10, trim = 0.025, epsilon = 1e-5, family = "gaussian") {
  n <- dim(data)[1]
  p <- dim(data)[2] - 1
  Fobj <- c(-beta, 0)
  cp_set <- list(NULL, 0)
  set <- c(0, 1)

  # choose the initial values based on pre-segmentation

  index <- rep(1:B, rep(n / B, B))
  coef.int <- matrix(NA, B, p)
  err_sd <- act_num <- rep(NA, B)
  for (i in 1:B)
  {
    cvfit <- cv.glmnet(as.matrix(data[index == i, 1:p]), data[index == i, p + 1], family = family)
    coef.int[i, ] <- coef(cvfit, s = "lambda.1se")[-1]
    resi <- data[index == i, p + 1] - as.matrix(data[index == i, 1:p]) %*% as.numeric(coef.int[i, ])
    err_sd[i] <- sqrt(mean(resi^2))
    act_num[i] <- sum(abs(coef.int[i, ]) > 0)
  }
  err_sd_mean <- mean(err_sd) # only works if error sd is unchanged.
  act_num_mean <- mean(act_num)
  beta <- (act_num_mean + 1) * beta # seems to work but there might be better choices

  X1 <- data[1, 1:p]
  cum_coef <- coef <- matrix(coef.int[1, ], p, 1)
  eta <- coef %*% X1
  # c_int <- diag(1/epsilon,p) - X1%o%X1/epsilon^2/(1+sum(X1^2)/epsilon)
  # cmatrix_inv <- array(c_int, c(p,p,1))
  cmatrix <- array(X1 %o% X1 + epsilon * diag(1, p), c(p, p, 1))

  for (t in 2:n)
  {
    m <- length(set)
    cval <- rep(NA, m)

    for (i in 1:(m - 1))
    {
      coef_c <- coef[, i]
      cum_coef_c <- cum_coef[, i]
      # cmatrix_inv_c <- cmatrix_inv[,,i]
      cmatrix_c <- cmatrix[, , i]
      k <- set[i] + 1
      out <- cost_lasso_update(data[t, ], coef_c, cum_coef_c, cmatrix_c, lambda = err_sd_mean * sqrt(2 * log(p) / (t - k + 1)))
      coef[, i] <- out[[1]]
      cum_coef[, i] <- out[[2]]
      # cmatrix_inv[,,i] <- out[[3]]
      cmatrix[, , i] <- out[[3]]
      if (t - k >= 2) cval[i] <- neg_log_lik_lasso(data[k:t, ], cum_coef[, i] / (t - k + 1), lambda = err_sd_mean * sqrt(2 * log(p) / (t - k + 1))) else cval[i] <- 0
    }

    # the choice of initial values requires further investigation

    cval[m] <- 0
    Xt <- data[t, 1:p]
    cum_coef_add <- coef_add <- coef.int[index[t], ]
    # cmatrix_inv_add <- diag(1/epsilon,p) - Xt%o%Xt/epsilon^2/(1+sum(Xt^2)/epsilon)

    coef <- cbind(coef, coef_add)
    cum_coef <- cbind(cum_coef, cum_coef_add)
    # cmatrix_inv <- abind::abind(cmatrix_inv, cmatrix_inv_add, along=3)
    cmatrix <- abind::abind(cmatrix, Xt %o% Xt + epsilon * diag(1, p), along = 3)

    obj <- cval + Fobj[set + 1] + beta
    min_val <- min(obj)
    ind <- which(obj == min_val)[1]
    cp_set_add <- c(cp_set[[set[ind] + 1]], set[ind])
    cp_set <- append(cp_set, list(cp_set_add))
    ind2 <- (cval + Fobj[set + 1]) <= min_val
    set <- c(set[ind2], t)
    coef <- coef[, ind2, drop = FALSE]
    cum_coef <- cum_coef[, ind2, drop = FALSE]
    cmatrix <- cmatrix[, , ind2, drop = FALSE]
    # cmatrix_inv <- cmatrix_inv[,,ind2,drop=FALSE]
    Fobj <- c(Fobj, min_val)
  }

  # Remove change-points close to the boundaries and merge change-points

  cp <- cp_set[[n + 1]]
  if (length(cp) > 0) {
    ind3 <- (1:length(cp))[(cp < trim * n) | (cp > (1 - trim) * n)]
    if (length(ind3) > 0) cp <- cp[-ind3]
  }

  cp <- sort(unique(c(0, cp)))
  index <- which((diff(cp) < trim * n) == TRUE)
  if (length(index) > 0) cp <- floor((cp[-(index + 1)] + cp[-index]) / 2)
  cp <- cp[cp > 0]

  # nLL <- 0
  # cp_loc <- unique(c(0,cp,n))
  # for(i in 1:(length(cp_loc)-1))
  # {
  #  seg <- (cp_loc[i]+1):cp_loc[i+1]
  #  data_seg <- data[seg,]
  #  out <- fastglm(as.matrix(data_seg[, 1:p]), data_seg[, p+1], family="binomial")
  #  nLL <- out$deviance/2 + nLL
  # }

  output <- list(cp)
  names(output) <- c("cp")
  return(output)
}

# Generate data from penalized linear regression models with change-points
#' @param n Number of observations.
#' @param d Dimension of the covariates.
#' @param true.coef True regression coefficients.
#' @param true.cp.loc True change-point locations.
#' @param Sigma Covariance matrix of the covariates.
#' @param evar Error variance.
#' @keywords internal
#'
#' @noRd
#' @return A list containing the generated data and the true cluster
#'    assignments.
data_gen_lasso <- function(n, d, true.coef, true.cp.loc, Sigma, evar) {
  loc <- unique(c(0, true.cp.loc, n))
  if (dim(true.coef)[2] != length(loc) - 1) stop("true.coef and true.cp.loc do not match")
  x <- mvtnorm::rmvnorm(n, mean = rep(0, d), sigma = Sigma)
  y <- NULL
  for (i in 1:(length(loc) - 1))
  {
    Xb <- x[(loc[i] + 1):loc[i + 1], , drop = FALSE] %*% true.coef[, i, drop = FALSE]
    add <- Xb + rnorm(length(Xb), sd = sqrt(evar))
    y <- c(y, add)
  }
  data <- cbind(x, y)
  true_cluster <- rep(1:(length(loc) - 1), diff(loc))
  result <- list(data, true_cluster)
  return(result)
}

testthat::test_that("penalized linear regression", {
  testthat::skip("This test is intended to be run manually.")
  set.seed(1)
  n <- 1000
  s <- 3
  d <- 50
  evar <- 0.5
  Sigma <- diag(1, d)
  true.cp.loc <- c(100, 300, 500, 800, 900)
  seg <- length(true.cp.loc) + 1
  true.coef <- matrix(rnorm(seg * s), s, seg)
  true.coef <- rbind(true.coef, matrix(0, d - s, seg))
  out <- data_gen_lasso(n, d, true.coef, true.cp.loc, Sigma, evar)
  data <- out[[1]]
  beta <- log(n) / 2 # beta here has different meaning

  change_points_lasso_fastcpd <- fastcpd::fastcpd(
    formula = y ~ . - 1,
    data = data.frame(y = data[, d + 1], x = data[, 1:d]),
    family = "lasso"
  )@cp_set

  change_points_lasso_fastcpd_sanity <- segd_lasso(data, beta, B = 10, trim = 0.025)$cp

  # TODO(doccstat): Deal with the bugs.
  # change_points_lasso_fastcpd_vanilla <- fastcpd::fastcpd(
  #   formula = y ~ . - 1,
  #   data = data.frame(y = data[, d + 1], x = data[, 1:d]),
  #   family = "lasso",
  #   vanilla_percentage = 1
  # )@cp_set

  change_points_lasso_fastcpd_vanilla_sanity <- pelt_vanilla_lasso(data, beta, cost = cost_lasso)$cp


  testthat::expect_equal(
    change_points_lasso_fastcpd,
    c(315, 798)
  )

  # TODO(doccstat): Check why the results are different.
  testthat::expect_equal(
    change_points_lasso_fastcpd_sanity,
    c(100, 300, 520, 800, 901)
  )

  # testthat::expect_equal(
  #   change_points_lasso_fastcpd_vanilla,
  #   c(374, 752, 1133)
  # )

  testthat::expect_equal(
    c(0, 103, 299, 510, 800, 900),
    change_points_lasso_fastcpd_vanilla_sanity
  )
})

# nolint end

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
