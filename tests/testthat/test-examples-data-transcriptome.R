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
        178, 265, 396, 529, 541, 578, 657, 789, 811, 870,
        935, 960, 1051, 1141, 1286, 1320, 1387, 1569, 1657,
        1725, 1907, 1973, 1994, 2041, 2059, 2142, 2201
      )
    )
  }
)
