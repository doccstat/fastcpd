testthat::test_that(
  "examples/fastcpd_meanvariance_2.R", {
    testthat::skip_if_not_installed("mvtnorm")
    source("examples/fastcpd_meanvariance_2.R")
    testthat::expect_equal(
      result@cp_set, c(2e+5, 3e+5, 5e+5, 7e+5, 8e+5), tolerance = 2
    )
  }
)
