testthat::test_that(
  "examples/data-well_log_2.R", {
    testthat::skip_if_not_installed("ggplot2")

    source("examples/data-well_log_2.R")
    testthat::expect_equal(
      result@cp_set,
      c(1021, 1057, 1502, 1661, 1842, 2023, 2385, 2445, 2507, 2567, 2744)
    )
  }
)
