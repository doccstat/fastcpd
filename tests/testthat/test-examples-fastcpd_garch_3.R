testthat::test_that(
  "examples/fastcpd_garch_3.txt", {
    testthat::skip_if_not_installed("ggplot2")

    examples_garch_3 <- readLines("examples/fastcpd_garch_3.txt")
    source(textConnection(paste(
      examples_garch_3[seq_len(length(examples_garch_3) - 2) + 1],
      collapse = "\n"
    )))

    testthat::expect_equal(result@cp_set, 74)
    testthat::expect_equal(max(result@residuals, na.rm = TRUE), 2.784562259)
    testthat::expect_equal(result@thetas, data.frame(
      "segment 1" = c(2.841266842, 0.135817457, 0.748930388),
      "segment 2" = c(0.180061756, 0.005695202, 0),
      check.names = FALSE
    ))
  }
)
