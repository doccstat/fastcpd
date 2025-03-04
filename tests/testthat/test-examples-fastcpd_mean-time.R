testthat::test_that(
  "examples/fastcpd_mean-time.R", {
    testthat::skip_on_cran()

    source("examples/fastcpd_mean-time.R")
    testthat::expect_equal(result@cp_set, 5e+6)
    testthat::expect_lt(result_time[3], 20)
  }
)
