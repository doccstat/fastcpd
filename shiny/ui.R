library(shiny)

ui <- fluidPage(

  theme = bslib::bs_theme(bootswatch = "flatly"),

  # App title
  titlePanel(
    paste(
      "fastcpd - Fast Change Point Detection",
      "(for demonstration purposes and reasonable data sizes only)"
    )
  ),

  selectInput(
    "family",
    "Choose a model family:",
    list(
      `Unlabeled Data` = list(
        "Mean Change" = "mean",
        "Variance Change" = "variance",
        "Mean and/or Variance Change" = "meanvariance"
      ),
      `Time Series` = list(
        "AR(p)" = "ar",
        "ARMA(p, q)" = "arma",
        "ARIMA(p, d, q)" = "arima",
        "GARCH(p, q)" = "garch",
        "VAR(p)" = "var"
      ),
      `Regression Data` = list(
        "Linear Regression" = "lm",
        "Penalized Linear Regression" = "lasso",
        "Logistic Regression" = "binomial",
        "Poisson Regression" = "poisson"
      ),
      `DIY` = list(
        "Custom" = "custom"
      )
    )
  ),

  conditionalPanel(
    condition = "input.family == 'custom'",
    textInput(
      inputId = "formula",
      label = "Enter a formula (e.g., y ~ . - 1):",
      value = "y ~ . - 1"
    )
  ),

  conditionalPanel(
    condition = "input.family == 'custom'",
    textInput(
      inputId = "custom_family",
      label = paste(
        "Enter a family (e.g., mean).",
        "Please keep the 'custom' family",
        "if cost functions will be provided below:"
      ),
      value = "custom"
    )
  ),

  conditionalPanel(
    condition = "input.family == 'custom'",
    textAreaInput(
      inputId = "cost",
      label = paste(
        "Enter the cost function",
        "(e.g., function(data) wiht optional `theta` parameter):"
      ),
      value = r"[
function(data, theta) {
  x <- data[, -1]
  y <- data[, 1]
  u <- x %*% theta
  nll <- -y * u + log(1 + exp(u))
  nll[u > 10] <- -y[u > 10] * u[u > 10] + u[u > 10]
  sum(nll)
}]",
      width = "600px",
      height = "250px"
    )
  ),

  conditionalPanel(
    condition = "input.family == 'custom'",
    textAreaInput(
      inputId = "cost_gradient",
      label = paste(
        "Enter the cost gradient function",
        "(e.g., function(data, theta)):"
      ),
      value = r"[
function(data, theta) {
  x <- data[nrow(data), -1]
  y <- data[nrow(data), 1]
  c(-(y - 1 / (1 + exp(-x %*% theta)))) * x
}]",
      width = "600px",
      height = "250px"
    )
  ),

  conditionalPanel(
    condition = "input.family == 'custom'",
    textAreaInput(
      inputId = "cost_hessian",
      label = paste(
        "Enter the cost Hessian function",
        "(e.g., function(data, theta)):"
      ),
      value = r"[
function(data, theta) {
  x <- data[nrow(data), -1]
  prob <- 1 / (1 + exp(-x %*% theta))
  (x %o% x) * c((1 - prob) * prob)
}]",
      width = "600px",
      height = "250px"
    )
  ),

  textInput(
    inputId = "beta",
    label = paste0(
      "Enter penalty selection criterion (BIC, MBIC or MDL) or a numeric value:"
    ),
    value = "MBIC"
  ),

  # Sidebar layout for file upload
  sidebarLayout(
    sidebarPanel(
      fileInput(
        inputId = "file",
        label = "Select a csv file (with header)",
        buttonLabel = "Select",
        multiple = FALSE,
        accept = c(".csv")
      )
    ),

    mainPanel(verbatimTextOutput(outputId = "summary"))
  )
)
