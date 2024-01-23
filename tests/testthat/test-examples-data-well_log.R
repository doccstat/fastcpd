testthat::test_that(
  "examples/data-well_log.R", {
    testthat::skip_if_not_installed("ggplot2")

    source("examples/data-well_log.R")
    testthat::expect_equal(
      result@cp_set,
      c(
        10, 572, 704, 830, 1021, 1057, 1198, 1347,
        1406, 1502, 1660, 1842, 2023, 2385, 2445,
        2507, 2567, 2749, 3709, 3820, 3976
      )
    )
  }
)
