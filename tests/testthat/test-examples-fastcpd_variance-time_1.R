testthat::test_that(
  "examples/fastcpd_variance-time_1.txt", {
    testthat::skip_on_cran()

    examples_variance_time <- readLines("examples/fastcpd_variance-time_1.txt")
    source(textConnection(paste(
      examples_variance_time[seq_len(length(examples_variance_time) - 2) + 1],
      collapse = "\n"
    )))

    testthat::expect_equal(result@cp_set, c(2997, 5989))
    testthat::expect_lt(result_time[3], 50)
  }
)
