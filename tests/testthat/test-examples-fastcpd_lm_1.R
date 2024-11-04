testthat::test_that(
  "examples/fastcpd_lm_1.R", {
    testthat::skip_if_not_installed("ggplot2")
    testthat::skip_if_not_installed("mvtnorm")

    source("examples/fastcpd_lm_1.R")
    testthat::expect_equal(result_lm@cp_set, c(97, 201))
  }
)
