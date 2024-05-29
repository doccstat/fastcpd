testthat::test_that(
  "examples/fastcpd_ar_1.R", {
    testthat::skip_if_not_installed("ggplot2")

    source("examples/fastcpd_ar_1.R")
    testthat::expect_equal(result@cp_set, 614)
  }
)
