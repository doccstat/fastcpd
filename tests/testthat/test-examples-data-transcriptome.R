testthat::test_that(
  "examples/data-transcriptome.txt", {
    examples_transcriptome <- readLines("examples/data-transcriptome.txt")
    source(textConnection(paste(
      examples_transcriptome[seq_len(length(examples_transcriptome) - 2) + 1],
      collapse = "\n"
    )))

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
