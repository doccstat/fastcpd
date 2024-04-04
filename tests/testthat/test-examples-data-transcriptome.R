testthat::test_that(
  "examples/data-transcriptome.txt", {
    testthat::skip_if_not_installed("ggplot2")
    testthat::skip_if_not_installed("gridExtra")

    examples_transcriptome <- readLines("examples/data-transcriptome.txt")
    source(textConnection(paste(
      examples_transcriptome[seq_len(length(examples_transcriptome) - 2) + 1],
      collapse = "\n"
    )))

    testthat::expect_equal(
      result@cp_set,
      c(
        177, 264, 401, 534, 578, 656, 788, 811, 869, 934, 960, 1051, 1141, 1286,
        1319, 1367, 1566, 1657, 1724, 1906, 1972, 1996, 2041, 2143, 2200
      )
    )
  }
)
