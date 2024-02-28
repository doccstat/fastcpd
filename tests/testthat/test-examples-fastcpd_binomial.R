testthat::test_that(
  "examples/fastcpd_binomial.R", {
    testthat::skip_if_not_installed("ggplot2")

    source("examples/fastcpd_binomial.R")
    testthat::expect_equal(result@cp_set, 301)
  }
)
