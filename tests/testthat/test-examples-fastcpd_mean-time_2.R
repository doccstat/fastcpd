testthat::test_that(
  "examples/fastcpd_mean-time_2.txt", {
    testthat::skip_on_cran()

    examples_mean_time <- readLines("examples/fastcpd_mean-time_2.txt")
    source(textConnection(paste(
      examples_mean_time[seq_len(length(examples_mean_time) - 2) + 1],
      collapse = "\n"
    )))

    testthat::expect_equal(result@cp_set, 5e+6)
    testthat::expect_lt(result_time[3], 20)
  }
)
