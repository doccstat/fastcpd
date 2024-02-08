library(shiny)

# Define UI for the app
ui <- fluidPage(

  theme = bslib::bs_theme(bootswatch = "flatly"),

  # App title
  titlePanel(
    paste0(
      "fastcpd - Fast Change Point Detection ",
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
        "MA(p)" = "ma",
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
      label = "Enter a family (e.g., mean):",
      value = "custom"
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

server <- function(input, output) {

  output$summary <- renderPrint(
    {
      input_file <- input$file
      if (!is.null(input_file)) {
        data <- read.csv(file = input_file$datapath)

        # Extract penalty selection criterion (beta) or a numeric value
        beta <- ifelse(
          input$beta %in% c("BIC", "MBIC", "MDL"),
          input$beta,
          suppressWarnings(as.numeric(input$beta))
        )

        # Use `fastcpd` wrapper functions so that a formula is not required.
        if (input$family == "custom") {
          result <- fastcpd::fastcpd(
            formula = eval(parse(text = input$formula)),
            data = data,
            beta = beta,
            family = input$custom_family
          )
        } else {
          fastcpd_wrapper <- switch(
            input$family,
            "mean" = fastcpd::fastcpd_mean,
            "variance" = fastcpd::fastcpd_variance,
            "meanvariance" = fastcpd::fastcpd_meanvariance,
            "ar" = fastcpd::fastcpd_ar,
            "ma" = fastcpd::fastcpd_ma,
            "arma" = fastcpd::fastcpd_arma,
            "arima" = fastcpd::fastcpd_arima,
            "garch" = fastcpd::fastcpd_garch,
            "var" = fastcpd::fastcpd_var,
            "lm" = fastcpd::fastcpd_lm,
            "lasso" = fastcpd::fastcpd_lasso,
            "binomial" = fastcpd::fastcpd_binomial,
            "poisson" = fastcpd::fastcpd_poisson
          )

          # Run the fastcpd wrapper function
          result <- fastcpd_wrapper(
            data = data,
            beta = beta
          )
        }

        summary(result)
      }
    },
  )

}

shinyApp(ui = ui, server = server)
