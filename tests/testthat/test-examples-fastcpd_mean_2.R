testthat::test_that(
  "examples/fastcpd_mean_2.R", {
    source("examples/fastcpd_mean_2.R")
    testthat::expect_equal(result@cp_set, c(3e+5, 7e+5))
  }
)
