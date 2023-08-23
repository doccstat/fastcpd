test_that("logistic regression", {
  # This is the same example with `fastcpd` documentation. Please keep it in
  # sync if the documentation ever changes.
  set.seed(1)

  kChangePointLocation <- 125
  kNumberOfDataPoints <- 300
  kDimension <- 5

  # There are kNumberOfDataPoints five-dimensional data points.
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

  change_points_binomial_fastcpd_sanity <- segd_binomial(
    cbind(x, y), (kDimension + 1) * log(kNumberOfDataPoints) / 2, B = 5
  )$cp
  expect_equal(change_points_binomial_fastcpd_sanity, kChangePointLocation)

  change_points_binomial_fastcpd_vanilla_sanity <- pelt_vanilla_binomial(
    cbind(x, y), (kDimension + 1) * log(kNumberOfDataPoints) / 2
  )$cp
  expect_equal(change_points_binomial_fastcpd_vanilla_sanity, c(0, kChangePointLocation))
})
