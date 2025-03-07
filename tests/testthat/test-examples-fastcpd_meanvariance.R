testthat::test_that(
  "examples/fastcpd_meanvariance.R", {
    testthat::skip_if_not_installed("mvtnorm")

    source("examples/fastcpd_meanvariance.R")
    testthat::expect_equal(result@cp_set, c(2e+5, 3e+5, 500408, 7e+5, 8e+5))
    testthat::expect_lt(result_time[3], 20)
  }
)
