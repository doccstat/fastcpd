testthat::test_that(
  "examples/fastcpd_garch_2.txt", {
    testthat::skip_if_not_installed("ggplot2")

    examples_garch_2 <- readLines("examples/fastcpd_garch_2.txt")
    source(textConnection(paste(
      examples_garch_2[seq_len(length(examples_garch_2) - 2) + 1],
      collapse = "\n"
    )))

    testthat::expect_equal(result@cp_set, 60)
    testthat::expect_equal(max(result@residuals, na.rm = TRUE), 2.567000894)
    testthat::expect_equal(result@thetas, data.frame(
      "segment 1" = c(13.834633386, 0.227888554, 0.020084797),
      "segment 2" = c(0.174437968, 0, 0),
      check.names = FALSE
    ))
  }
)
