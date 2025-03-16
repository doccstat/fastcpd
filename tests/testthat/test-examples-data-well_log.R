testthat::test_that(
  "examples/data-well_log.R", {
    testthat::skip_if_not_installed("ggplot2")

    source("examples/data-well_log.R")
    testthat::expect_equal(
      result@cp_set,
      c(
        7, 19, 356, 448, 717, 844, 1034, 1070, 1215, 1369, 1428,
        1526, 1685, 1866, 2047, 2409, 2469, 2531, 2591, 2777,
        3490, 3533, 3672, 3744, 3855, 3886, 3945, 3963, 4035
      )
    )
  }
)
