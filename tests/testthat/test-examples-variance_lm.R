testthat::test_that(
  "examples/variance_lm.R", {
    testthat::skip_if_not_installed("mvtnorm")

    source("examples/variance_lm.R")
    testthat::expect_equal(sigma2, 10.4922713)
    testthat::expect_equal(
      sigma, matrix(c(3.71946553, 0.03459023, 0.01962422, 3.59825273), 2, 2)
    )
  }
)
