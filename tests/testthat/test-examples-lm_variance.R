testthat::test_that(
  "examples/lm_variance.R", {
    testthat::skip_if_not_installed("mvtnorm")

    source("examples/lm_variance.R")
    testthat::expect_equal(sigma2, 10.4922713)
  }
)
