testthat::test_that(
  "examples/fastcpd_ts.txt", {
    testthat::skip_if_not_installed("ggplot2")

    examples_ts <- readLines("examples/fastcpd_ts.txt")
    source(textConnection(paste(
      examples_ts[seq_len(length(examples_ts) - 2) + 1],
      collapse = "\n"
    )))

    testthat::expect_equal(result@cp_set, c(203, 356))
  }
)
