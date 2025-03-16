testthat::test_that(
  "examples/fastcpd_ar_1.R", {
    testthat::skip_if_not_installed("ggplot2")

    source("examples/fastcpd_ar_1.R")
    testthat::expect_equal(result@cp_set, 614)
    testthat::expect_equal(result@thetas, data.frame(
      "segment 1" = c(0.57120256, -0.20985108, 0.08221978),
      "segment 2" = c(0.237180884, 0.403124372, 0.229032322),
      check.names = FALSE
    ))
  }
)
