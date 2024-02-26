testthat::test_that(
  "examples/variance_lm.R", {
    testthat::skip_if_not_installed("mvtnorm")

    source("examples/variance_lm.R")
    testthat::expect_equal(sigma2, 10.4922713)
  }
)
