testthat::test_that(
  "examples/fastcpd_mean-time_1.R", {
    testthat::skip_on_cran()

    source("examples/fastcpd_mean-time_1.R")
    testthat::expect_equal(result@cp_set, 10007)
    testthat::expect_lt(result_time[3], 30)
  }
)
