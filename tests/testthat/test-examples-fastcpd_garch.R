testthat::test_that(
  "examples/fastcpd_garch.txt", {
    testthat::skip_if_not_installed("ggplot2")

    examples_garch <- readLines("examples/fastcpd_garch.txt")
    source(textConnection(paste(
      examples_garch[seq_len(length(examples_garch) - 2) + 1],
      collapse = "\n"
    )))

    testthat::expect_equal(result@cp_set, 759)
    testthat::expect_equal(max(result@residuals, na.rm = TRUE), 3.89822750)
    testthat::expect_equal(result@thetas, data.frame(
      "segment 1" = c(79.822854422, 0.196544235, 0.000000000),
      "segment 2" = c(1.376189971, 0.125421192, 0.356561246),
      check.names = FALSE
    ))
  }
)
