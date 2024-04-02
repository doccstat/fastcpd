testthat::test_that(
  "examples/fastcpd_mean-time_3.R", {
    source("examples/fastcpd_mean-time_3.R")
    testthat::expect_equal(result@cp_set, c(10034, 20082))
    testthat::expect_lt(result_time[3], 100)
  }
)
