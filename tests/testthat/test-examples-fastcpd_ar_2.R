testthat::test_that(
  "examples/fastcpd_ar_2.R", {
    testthat::skip_if_not_installed("ggplot2")

    source("examples/fastcpd_ar_2.R")
    testthat::expect_equal(result@cp_set, 614)
    testthat::expect_lt(
      abs(median(result@residuals, na.rm = TRUE) + 0.11998684), 1e-8
    )
  }
)
