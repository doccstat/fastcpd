testthat::test_that(
  "examples/data-well_log_1.R", {
    testthat::skip_if_not_installed("ggplot2")

    source("examples/data-well_log_1.R")
    testthat::expect_equal(
      result@cp_set,
      c(
        12, 566, 704, 776, 1021, 1057, 1198, 1347, 1406, 1502,
        1665, 1842, 2023, 2202, 2384, 2445, 2507, 2567, 2749,
        2926, 3094, 3107, 3509, 3622, 3709, 3820, 3976
      )
    )
  }
)
