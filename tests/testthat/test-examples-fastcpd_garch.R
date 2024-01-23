testthat::test_that(
  "examples/fastcpd_garch.txt", {
    testthat::skip_if_not_installed("ggplot2")

    examples_garch <- readLines("examples/fastcpd_garch.txt")
    source(textConnection(paste(
      examples_garch[seq_len(length(examples_garch) - 2) + 1],
      collapse = "\n"
    )))

    testthat::expect_equal(result@cp_set, 205)
  }
)
