testthat::test_that(
  "examples/fastcpd_garch_2.R", {
    testthat::skip_if_not_installed("ggplot2")

    source("examples/fastcpd_garch_2.R")

    testthat::expect_equal(result@cp_set, 98)
    testthat::expect_equal(max(result@residuals, na.rm = TRUE), 2.603021577)
    testthat::expect_equal(result@thetas, data.frame(
      "segment 1" = c(5.671883227, 0.183053584, 0.565126266),
      "segment 2" = c(0.194429864, 0.081632616, 0.041330328),
      check.names = FALSE
    ))
  }
)
