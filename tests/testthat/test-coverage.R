testthat::test_that(
  "1d mean", {
    set.seed(1)
    data <- c(rnorm(300, 0, 100), rnorm(400, 100, 100), rnorm(300, 0, 100))
    result_mean <- fastcpd.mean(data)
    testthat::expect_equal(result_mean@cp_set, c(294, 702))
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
      cost = logistic_loss,
      cost_gradient = logistic_loss_gradient,
      cost_hessian = logistic_loss_hessian,
      k = function(x) if (x < 10) 1 else 0
    )
    testthat::expect_equal(result@cp_set, 147)
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
        data.frame(y = y, x = x), vanilla_percentage = 1, warm_start = TRUE
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

testthat::test_that("beta x2", {
  testthat::expect_equal(
    fastcpd.mean(well_log, beta = log(length(well_log)))@cp_set,
    c(566, 740, 1039, 1198, 1424, 1661, 1842, 2023, 2476, 2744, 3709, 3820)
  )
})
