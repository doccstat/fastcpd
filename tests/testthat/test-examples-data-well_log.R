testthat::test_that(
  "examples/data-well_log.R", {
    testthat::skip_if_not_installed("ggplot2")

    source("examples/data-well_log.R")
    testthat::expect_equal(
      result@cp_set,
      c(
        12, 1031, 1056, 1198, 1407, 1504, 1661,
        1842, 2024, 2385, 2445, 2508, 2567, 2747
      )
    )
  }
)
