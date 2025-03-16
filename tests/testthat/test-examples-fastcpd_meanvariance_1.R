testthat::test_that(
  "examples/fastcpd_meanvariance_1.R", {
    source("examples/fastcpd_meanvariance_1.R")
    testthat::expect_equal(result@cp_set, c(3000, 4000, 7000))
  }
)
