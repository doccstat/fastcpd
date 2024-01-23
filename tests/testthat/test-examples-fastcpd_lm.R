testthat::test_that(
  "examples/fastcpd_lm.R", {
    testthat::skip_if_not_installed("ggplot2")
    testthat::skip_if_not_installed("mvtnorm")

    source("examples/fastcpd_lm.R")
    testthat::expect_equal(result@cp_set, c(97, 201))
  }
)
