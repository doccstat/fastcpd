testthat::test_that(
  "examples/fastcpd_meanvariance_1.R", {
    testthat::skip_if_not_installed("mvtnorm")

    source("examples/fastcpd_meanvariance_1.R")
    testthat::expect_equal(result@cp_set, c(300, 700, 1001, 1300, 1700))
  }
)
