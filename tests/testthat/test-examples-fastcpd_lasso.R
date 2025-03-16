testthat::test_that(
  "examples/fastcpd_lasso.txt", {
    testthat::skip_if_not_installed("ggplot2")
    testthat::skip_if_not_installed("mvtnorm")

    examples_lasso <- readLines("examples/fastcpd_lasso.txt")
    source(textConnection(paste(
      examples_lasso[seq_len(length(examples_lasso) - 2) + 1],
      collapse = "\n"
    )))

    testthat::expect_equal(result@cp_set, c(80, 200, 321))
  }
)
