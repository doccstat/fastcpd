testthat::test_that(
  "examples/fastcpd_4.txt", {
    testthat::skip_if_not_installed("ggplot2")
    testthat::skip_if_not_installed("mvtnorm")

    source("examples/fastcpd_4.R")
    testthat::expect_equal(huber_regression_result@cp_set, c(401, 726))
  }
)
