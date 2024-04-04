library(shiny)

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
          cost <- NULL
          cost_gradient <- NULL
          cost_hessian <- NULL
          if (input$custom_family == "custom") {
            cost <- eval(parse(text = input$cost))
            cost_gradient <- eval(parse(text = input$cost_gradient))
            cost_hessian <- eval(parse(text = input$cost_hessian))
          }
          result <- fastcpd::fastcpd(
            formula = eval(parse(text = input$formula)),
            data = data,
            beta = beta,
            family = input$custom_family,
            cost = cost,
            cost_gradient = cost_gradient,
            cost_hessian = cost_hessian
          )
        } else {
          fastcpd_wrapper <- switch(
            input$family,
            "mean" = fastcpd::fastcpd_mean,
            "variance" = fastcpd::fastcpd_variance,
            "meanvariance" = fastcpd::fastcpd_meanvariance,
            "ar" = fastcpd::fastcpd_ar,
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
