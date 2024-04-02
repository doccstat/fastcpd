testthat::test_that(
  "examples/fastcpd_mean-time_3.txt", {
    examples_mean_time <- readLines("examples/fastcpd_mean-time_3.txt")
    source(textConnection(paste(
      examples_mean_time[seq_len(length(examples_mean_time) - 2) + 1],
      collapse = "\n"
    )))

    testthat::expect_equal(result@cp_set, c(10034, 20082))
    testthat::expect_lt(result_time[3], 100)
  }
)
