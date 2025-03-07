testthat::test_that(
  "examples/fastcpd_variance_1.R", {
    source("examples/fastcpd_variance_1.R")
    testthat::expect_equal(result@cp_set, c(3e+5, 7e+5))
    testthat::expect_lt(result_time[3], 50)
  }
)
