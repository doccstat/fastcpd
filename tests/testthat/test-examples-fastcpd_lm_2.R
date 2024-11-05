testthat::test_that(
  "examples/fastcpd_lm_2.txt", {
    testthat::skip_if_not_installed("mvtnorm")

    examples_lm <- readLines("examples/fastcpd_lm_2.txt")
    source(textConnection(paste(
      examples_lm[seq_len(length(examples_lm) - 2) + 1],
      collapse = "\n"
    )))

    testthat::expect_equal(result_mlm@cp_set, 350)
  }
)
