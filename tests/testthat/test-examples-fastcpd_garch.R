testthat::test_that(
  "examples/fastcpd_garch.txt", {
    testthat::skip_if(Sys.getenv("R_COVR") == "true", "slow GARCH example: skip during coverage")
    testthat::skip_if_not_installed("ggplot2")

    examples_garch <- readLines("examples/fastcpd_garch.txt")
    source(textConnection(paste(
      examples_garch[seq_len(length(examples_garch) - 2) + 1],
      collapse = "\n"
    )))

    testthat::expect_equal(result@cp_set, 756)
    testthat::expect_equal(max(result@residuals, na.rm = TRUE), 3.496100993)
    testthat::expect_equal(result@thetas, data.frame(
      "segment 1" = c(19.036319299, 0.639478846, 0.216478285),
      "segment 2" = c(1.516258347, 0.132955389, 0.296962417),
      check.names = FALSE
    ))
  }
)
