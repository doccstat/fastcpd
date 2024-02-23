testthat::test_that(
  "examples/fastcpd_mean-time.R", {
    source("examples/fastcpd_mean-time.R")
    testthat::expect_equal(result@cp_set, 10006)
    testthat::expect_true(result_time[3] < 2)
  }
)
