testthat::test_that(
  "examples/fastcpd_mean-time.R", {
    source("examples/fastcpd_mean-time.R")
    testthat::expect_equal(result@cp_set, 10006)

    # Something is seriously wrong if the running time is greater than
    # 30 seconds.
    testthat::expect_lt(result_time[3], 30)
  }
)
