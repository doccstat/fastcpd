testthat::test_that(
  "examples/fastcpd_garch.txt", {
    testthat::skip_if_not_installed("ggplot2")

    examples_garch <- readLines("examples/fastcpd_garch.txt")
    source(textConnection(paste(
      examples_garch[seq_len(length(examples_garch) - 2) + 1],
      collapse = "\n"
    )))

    testthat::expect_equal(result@cp_set, 205)
    testthat::expect_equal(max(result@residuals, na.rm = TRUE), 2.579871666)
    testthat::expect_equal(result@thetas, data.frame(
      "segment 1" = c(11.536369486, 0.469538255, 0.342805214),
      "segment 2" = c(1.998385819, 0.012825797, 0.136819886),
      check.names = FALSE
    ))
  }
)
