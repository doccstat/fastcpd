testthat::test_that(
  "examples/data-well_log-quantile.R", {
    testthat::skip_if_not_installed("ggplot2")
    testthat::skip_if_not_installed("matrixStats")

    examples_transcriptome <- readLines("examples/data-well_log-quantile.txt")
    source(textConnection(paste(
      examples_transcriptome[seq_len(length(examples_transcriptome) - 2) + 1],
      collapse = "\n"
    )))

    testthat::expect_equal(
      result@cp_set,
      c(
        19, 577, 1034, 1070, 1216, 1361, 1428, 1526, 1685, 1866, 2047,
        2409, 2469, 2531, 2591, 2775, 3744, 3855, 3945, 3963
      )
    )
  }
)
