testthat::test_that(
  "examples/fastcpd_variance.R", {
    testthat::skip_if_not_installed("mvtnorm")

    source("examples/fastcpd_variance.R")
    testthat::expect_equal(result@cp_set, c(300, 700))
  }
)
