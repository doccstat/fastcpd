testthat::test_that(
  "man/examples/fastcpd_ts.txt", {
    testthat::skip_if_not_installed("forecast")

    examples_ts <- readLines("man/examples/fastcpd_ts.txt")
    source(textConnection(paste(
      examples_ts[seq_len(length(examples_ts) - 2) + 1],
      collapse = "\n"
    )))

    testthat::expect_equal(result_ts@cp_set, 178)
  }
)
