test_that("summary output", {
  match_call <- ""
  class(match_call) <- "language"

  class_instance <- methods::new(
    Class = "fastcpd",
    call = match_call,
    data = data.frame(matrix(NA, 0, 0)),
    family = "custom",
    cp_set = c(100, 200),
    cost_values = c(9.5, 10.5),
    residuals = 0,
    thetas = data.frame(matrix(NA, 0, 0)),
    cp_only = TRUE
  )

  expect_equal(
    capture.output(summary(class_instance)),
    c("", "Call:", "structure(\"\", class = \"language\")", "", "Change points:", "100 200 ")
  )
})
