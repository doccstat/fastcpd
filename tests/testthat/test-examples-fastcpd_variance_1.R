testthat::test_that(
  "examples/fastcpd_variance_1.R", {
    source("examples/fastcpd_variance_1.R")
    testthat::expect_equal(result@cp_set, c(300, 700))
  }
)
