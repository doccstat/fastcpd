library(shiny)

# Define UI for the app
ui <- fluidPage(

  # App title
  titlePanel("fastcpd - Fast Change Point Detection"),

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

        # Use `fastcpd` wrapper functions so that a formula is not required.
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
          "poisson" = fastcpd::fastcpd_poisson,
          "custom" = fastcpd::fastcpd
        )
        result <- fastcpd_wrapper(
          data = data
        )
        summary(result)
      }
    },
  )

}

shinyApp(ui = ui, server = server)
