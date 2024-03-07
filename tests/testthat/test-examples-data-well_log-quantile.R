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
        566, 1021, 1057, 1340, 1502, 1661, 1842, 2023,
        2385, 2445, 2507, 2567, 2744, 3709, 3820
      )
    )
  }
)
