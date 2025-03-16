library(testthat)
library(fastcpd)

if (identical(Sys.getenv("NOT_CRAN"), "true")) {
  test_check("fastcpd")
}
