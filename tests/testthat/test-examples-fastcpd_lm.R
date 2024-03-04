testthat::test_that(
  "examples/fastcpd_lm.R", {
    testthat::skip_if_not_installed("ggplot2")
    testthat::skip_if_not_installed("mvtnorm")

    source("examples/fastcpd_lm.R")
    testthat::expect_equal(result_lm@cp_set, c(97, 201))
    testthat::expect_equal(result_mlm@cp_set, 375)
  }
)
