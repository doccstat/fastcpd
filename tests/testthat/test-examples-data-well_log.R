testthat::test_that(
  "examples/data-well_log.txt", {
    examples_well_log <- readLines("examples/data-well_log.txt")
    source(textConnection(paste(
      examples_well_log[seq_len(length(examples_well_log) - 2) + 1],
      collapse = "\n"
    )))

    testthat::expect_equal(
      result@cp_set,
      c(
        10, 566, 704, 776, 1021, 1057, 1198, 1347, 1406, 1502, 1661,
        1842, 2023, 2385, 2445, 2507, 2567, 2744, 3709, 3820, 3976
      )
    )
  }
)
