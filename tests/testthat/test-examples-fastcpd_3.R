testthat::test_that(
  "examples/fastcpd_3.txt", {
    examples_garch <- readLines("examples/fastcpd_3.txt")
    source(textConnection(paste(
      examples_garch[seq_len(length(examples_garch) - 2) + 1],
      collapse = "\n"
    )))

    testthat::expect_equal(result_builtin@cp_set, 198)
    testthat::expect_equal(result_custom@cp_set, 204)
  }
)
