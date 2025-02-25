testthat::test_that(
  "examples/fastcpd_arma.R", {
    testthat::skip_if_not_installed("ggplot2")

    examples_arima <- readLines("examples/fastcpd_arma.txt")
    source(textConnection(paste(
      examples_arima[seq_len(length(examples_arima) - 2) + 1],
      collapse = "\n"
    )))

    testthat::expect_equal(result@cp_set, 142)
  }
)
