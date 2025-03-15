testthat::test_that(
  "examples/fastcpd_mean_1.R", {
    source("examples/fastcpd_mean_1.R")
    testthat::expect_equal(result@cp_set, c(300, 700))
  }
)
