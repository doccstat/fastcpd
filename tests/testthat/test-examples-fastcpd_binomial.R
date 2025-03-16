testthat::test_that(
  "examples/fastcpd_binomial.R", {
    testthat::skip_if_not_installed("ggplot2")
    testthat::skip_if_not_installed("mvtnorm")

    source("examples/fastcpd_binomial.R")
    testthat::expect_equal(result@cp_set, 302)
    testthat::expect_equal(result@thetas, data.frame(
      "segment 1" = c(-0.926018205, -1.603383466, 1.034333794, 0.365386960),
      "segment 2" = c(2.129496165, 2.758324727, 2.381801028, 0.726115162),
      check.names = FALSE
    ))
  }
)
