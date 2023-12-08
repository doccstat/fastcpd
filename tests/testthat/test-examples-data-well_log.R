testthat::test_that(
  "examples/data-well_log.R", {
    source("examples/data-well_log.R")
    testthat::expect_equal(
      result@cp_set,
      c(
        10, 566, 704, 776, 1021, 1057, 1198, 1347, 1406, 1502,
        1660, 1842, 2023, 2202, 2384, 2445, 2507, 2567, 2749,
        2921, 3072, 3509, 3622, 3709, 3820, 3976
      )
    )
  }
)
