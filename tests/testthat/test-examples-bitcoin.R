testthat::test_that(
  "examples/bitcoin.txt", {
    examples_bitcoin <- readLines("examples/bitcoin.txt")
    source(textConnection(paste(
      examples_bitcoin[seq_len(length(examples_bitcoin) - 2) + 1],
      collapse = "\n"
    )))

    testthat::expect_equal(result@cp_set, 113)
  }
)
