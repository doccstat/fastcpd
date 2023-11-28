testthat::test_that(
  "examples/fastcpd_arma.R", {
    examples_arima <- readLines("examples/fastcpd_arma.txt")
    source(textConnection(paste(
      examples_arima[seq_len(length(examples_arima) - 2) + 1],
      collapse = "\n"
    )))

    # TODO(doccstat): Deal with the randomness in the example.
    testthat::expect_equal(
      result@cp_set, c(344, 385)[1 + (Sys.info()["sysname"] != "Darwin")]
    )
  }
)
