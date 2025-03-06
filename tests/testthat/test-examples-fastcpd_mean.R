testthat::test_that(
  "examples/fastcpd_mean.R", {
    testthat::skip_if_not_installed("mvtnorm")

    source("examples/fastcpd_mean.R")
    testthat::expect_equal(result@cp_set, c(3e+5, 7e+5))
    testthat::expect_lt(result_time[3], 20)
  }
)
