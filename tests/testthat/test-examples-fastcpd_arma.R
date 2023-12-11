testthat::skip("These tests are intended to be run manually.")

testthat::test_that(
  "examples/fastcpd_arma.R", {
    examples_arima <- readLines("examples/fastcpd_arma.txt")
    source(textConnection(paste(
      examples_arima[seq_len(length(examples_arima) - 2) + 1],
      collapse = "\n"
    )))

    # TODO(doccstat): Deal with the randomness in the example.
    # Local mac mini with M1 gives change point at 344.
    testthat::expect_equal(result@cp_set, 385)
  }
)
