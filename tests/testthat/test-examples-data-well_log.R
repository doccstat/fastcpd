testthat::test_that(
  "examples/data-well_log.R", {
    testthat::skip_if_not_installed("ggplot2")

    source("examples/data-well_log.R")
    testthat::expect_equal(
      result@cp_set,
      c(
        7, 19, 65, 356, 445, 717, 792, 1034, 1070, 1215, 1368, 1428,
        1526, 1684, 1866, 2047, 2409, 2469, 2531, 2591, 2775, 3166,
        3314, 3490, 3533, 3673, 3744, 3855, 3886, 3945, 3963, 4035
      )
    )
  }
)
