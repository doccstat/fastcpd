testthat::test_that(
  "examples/fastcpd_variance.R", {
    source("examples/fastcpd_variance.R")
    testthat::expect_equal(result@cp_set, c(300, 700))
  }
)
