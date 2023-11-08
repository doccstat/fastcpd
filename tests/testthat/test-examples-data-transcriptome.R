testthat::test_that(
  "examples/data-transcriptome.R", {
    source("examples/data-transcriptome.R")
    testthat::expect_equal(
      result@cp_set,
      c(
        135, 177, 264, 394, 531, 579, 656, 788, 811, 869, 934,
        966, 1051, 1141, 1253, 1286, 1319, 1368, 1568, 1657,
        1674, 1724, 1906, 1972, 1994, 2041, 2058, 2146, 2200
      )
    )
  }
)
