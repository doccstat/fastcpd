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
  "ar mdl", {
    set.seed(1)
    n <- 1000
    x <- rep(0, n + 3)
    for (i in 1:600) {
      x[i + 3] <- 0.6 * x[i + 2] - 0.2 * x[i + 1] + 0.1 * x[i] + rnorm(1, 0, 3)
    }
    for (i in 601:1000) {
      x[i + 3] <- 0.3 * x[i + 2] + 0.4 * x[i + 1] + 0.2 * x[i] + rnorm(1, 0, 3)
    }
    result <- fastcpd.ar(x[3 + seq_len(n)], 3)
    testthat::expect_equal(result@cp_set, 614)
    testthat::expect_lt(
      abs(median(result@residuals, na.rm = TRUE) + 0.11998684), 1e-8
    )
  }
)

testthat::test_that(
  "1d variance", {
    set.seed(1)
    data <- c(rnorm(300, 0, 1), rnorm(400, 0, 100), rnorm(300, 0, 1))
    result_variance <- fastcpd.variance(data)
    testthat::expect_equal(result_variance@cp_set, c(300, 700))
  }
)

testthat::test_that(
  "1d mv", {
    set.seed(1)
    data <- c(
      rnorm(300, 0, 1), rnorm(400, 10, 1), rnorm(300, 0, 50),
      rnorm(300, 0, 1), rnorm(400, 10, 1), rnorm(300, 10, 50)
    )
    result_meanvariance <- fastcpd.meanvariance(data)
    testthat::expect_equal(
      result_meanvariance@cp_set, c(300, 700, 1000, 1300, 1700)
    )
  }
)

testthat::test_that(
  "1d custom multiple epochs", {
    set.seed(1)
    p <- 1
    x <- matrix(rnorm(300 * p, 0, 1), ncol = p)
    theta <- rbind(rnorm(p, 0, 1), rnorm(p, 4, 1))
    y <- c(
      rbinom(150, 1, 1 / (1 + exp(-x[1:150, ] * theta[1, ]))),
      rbinom(150, 1, 1 / (1 + exp(-x[151:300, ] * theta[2, ])))
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
    result <- fastcpd(
      formula = y ~ . - 1,
      data = data.frame(y = y, x = x),
      beta = "BIC",
      cost = logistic_loss,
      cost_gradient = logistic_loss_gradient,
      cost_hessian = logistic_loss_hessian,
      multiple_epochs = function(segment_length) {
        if (segment_length < 10) 1 else 0
      },
      r.progress = FALSE
    )
    testthat::expect_equal(result@cp_set, 137)
  }
)

testthat::test_that(
  "2d custom", {
    set.seed(1)
    p <- 2
    x <- matrix(rnorm(300 * p, 0, 1), ncol = p)
    theta <- rbind(rnorm(p, 0, 1), rnorm(p, 4, 1))
    y <- c(
      rbinom(150, 1, 1 / (1 + exp(-x[1:150, ] * theta[1, ]))),
      rbinom(150, 1, 1 / (1 + exp(-x[151:300, ] * theta[2, ])))
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
    result <- fastcpd(
      formula = y ~ . - 1,
      data = data.frame(y = y, x = x),
      cost = logistic_loss,
      cost_gradient = logistic_loss_gradient,
      cost_hessian = logistic_loss_hessian
    )
    testthat::expect_equal(result@cp_set, 153)
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
