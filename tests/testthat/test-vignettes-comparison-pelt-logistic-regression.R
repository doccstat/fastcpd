# Everything in this script is provided as is. The purpose of this script is to
# do a sanity check on the C++ implementation of `fastcpd`.

testthat::skip_if_not_installed("mvtnorm")

# nolint start: script provided as is

testthat::test_that("logistic regression", {
  set.seed(1)
  p <- 5
  x <- matrix(rnorm(300 * p, 0, 1), ncol = p)

  # Randomly generate coefficients with different means.
  theta <- rbind(rnorm(p, 0, 1), rnorm(p, 2, 1))

  # Randomly generate response variables based on the segmented data and
  # corresponding coefficients
  y <- c(
    rbinom(125, 1, 1 / (1 + exp(-x[1:125, ] %*% theta[1, ]))),
    rbinom(300 - 125, 1, 1 / (1 + exp(-x[(125 + 1):300, ] %*% theta[2, ])))
  )

  change_points_binomial_fastcpd <- suppressWarnings(fastcpd.binomial(
    cbind(y, x),
    segment_count = 5,
    beta = "BIC",
    cost_adjustment = "BIC"
  ))@cp_set

  testthat::expect_equal(change_points_binomial_fastcpd, 125)

  warning_messages <- testthat::capture_warnings(
    change_points_binomial_fastcpd_vanilla <- fastcpd.binomial(
      cbind(y, x),
      segment_count = 5,
      vanilla_percentage = 1,
      beta = "BIC",
      cost_adjustment = "BIC"
    )@cp_set
  )

  testthat::expect_length(warning_messages, 3831)

  testthat::expect_equal(change_points_binomial_fastcpd_vanilla, 125)
})

# nolint end
