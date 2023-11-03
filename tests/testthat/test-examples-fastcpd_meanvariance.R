testthat::test_that(
  "examples/fastcpd_meanvariance.txt", {
    examples_meanvariance <- readLines("examples/fastcpd_meanvariance.txt")
    source(textConnection(paste(
      examples_meanvariance[seq_len(length(examples_meanvariance) - 2) + 1],
      collapse = "\n"
    )))

    testthat::expect_equal(result@cp_set, c(300, 700, 1000, 1300, 1700))
  }
)
