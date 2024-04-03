testthat::test_that(
  "examples/data-well_log.R", {
    testthat::skip_if_not_installed("ggplot2")

    source("examples/data-well_log.R")
    testthat::expect_equal(
      result@cp_set,
      c(
        13, 560, 705, 777, 1022, 1057, 1199, 1347, 1406, 1503, 1661,
        1843, 2024, 2385, 2446, 2508, 2568, 2749, 3710, 3820, 3977
      )
    )
  }
)
