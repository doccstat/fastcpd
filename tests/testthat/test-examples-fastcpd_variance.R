testthat::test_that(
  "examples/fastcpd_variance.R", {
    testthat::skip_if_not_installed("mvtnorm")
    source("examples/fastcpd_variance.R")
    testthat::expect_equal(result@cp_set, c(300002, 700003))
    testthat::expect_lt(result_time[3], 20)
  }
)
