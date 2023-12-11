testthat::test_that(
  "examples/fastcpd_poisson.txt", {
    testthat::skip_if_not_installed("mvtnorm")

    examples_poisson <- readLines("examples/fastcpd_poisson.txt")
    source(textConnection(paste(
      examples_poisson[seq_len(length(examples_poisson) - 2) + 1],
      collapse = "\n"
    )))

    testthat::expect_equal(result@cp_set, c(300, 704, 1003))
  }
)
