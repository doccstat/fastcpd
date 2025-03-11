testthat::test_that(
  "examples/fastcpd_meanvariance.R", {
    testthat::skip_if_not_installed("mvtnorm")

    source("examples/fastcpd_meanvariance.R")
    testthat::expect_equal(
      result@cp_set, c(2e+5, 3e+5, 5e+5, 7e+5, 8e+5), tolerance = 2
    )
    testthat::expect_lt(result_time[3], 50)
  }
)
