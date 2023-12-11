testthat::skip("These tests are intended to be run manually.")

testthat::test_that(
  "examples/fastcpd_arima.txt", {
    testthat::skip_if_not_installed("forecast")

    examples_arima <- readLines("examples/fastcpd_arima.txt")
    source(textConnection(paste(
      examples_arima[seq_len(length(examples_arima) - 2) + 1],
      collapse = "\n"
    )))

    # TODO(doccstat): Deal with the randomness in the example.
    testthat::expect_equal(
      result@cp_set, c(270, 317)[1 + (Sys.info()["sysname"] == "Windows")]
    )
  }
)
