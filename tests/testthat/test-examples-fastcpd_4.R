testthat::test_that(
  "examples/fastcpd_4.txt", {
    testthat::skip_if_not_installed("mvtnorm")

    examples_garch <- readLines("examples/fastcpd_4.txt")
    source(textConnection(paste(
      examples_garch[seq_len(length(examples_garch) - 2) + 1],
      collapse = "\n"
    )))

    testthat::expect_equal(result_no_vp@cp_set, c(79, 202, 325))
    testthat::expect_equal(result_20_vp@cp_set, c(80, 202, 320))
  }
)
