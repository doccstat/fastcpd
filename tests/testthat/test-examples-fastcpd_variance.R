testthat::test_that(
  "examples/fastcpd_variance.txt", {
    examples_mean <- readLines("examples/fastcpd_variance.txt")
    source(textConnection(paste(
      examples_mean[seq_len(length(examples_mean) - 2) + 1],
      collapse = "\n"
    )))

    testthat::expect_equal(result@cp_set, c(300, 700))
  }
)
