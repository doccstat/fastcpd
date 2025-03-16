testthat::test_that(
  "examples/plot.R", {
    testthat::skip_if_not_installed("ggplot2")
    testthat::skip_if_not_installed("mvtnorm")

    source("examples/plot.R")
    testthat::expect_equal(result@cp_set, c(100, 201))
  }
)
