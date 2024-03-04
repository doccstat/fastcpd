testthat::test_that(
  "examples/fastcpd_3.txt", {
    testthat::skip_if_not_installed("ggplot2")
    testthat::skip_if_not_installed("mvtnorm")

    examples_garch <- readLines("examples/fastcpd_3.txt")
    source(textConnection(paste(
      examples_garch[seq_len(length(examples_garch) - 2) + 1],
      collapse = "\n"
    )))

    testthat::expect_equal(result_builtin@cp_set, 200)
    testthat::expect_equal(result_custom@cp_set, 201)
    testthat::expect_equal(result_custom_two_epochs@cp_set, 200)
  }
)
