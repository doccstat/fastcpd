testthat::test_that(
  "examples/fastcpd_ma.txt", {
    testthat::skip_if_not_installed("forecast")

    examples_ma <- readLines("examples/fastcpd_ma.txt")
    source(textConnection(paste(
      examples_ma[seq_len(length(examples_ma) - 2) + 1],
      collapse = "\n"
    )))

    testthat::expect_equal(result@cp_set, 200)
  }
)
