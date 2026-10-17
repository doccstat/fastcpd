# Detailed plots show the first response and its matching residual column.
plot_x <- seq_len(16)
plot_response <- 2 * plot_x + sin(plot_x)
plot_residual_fit <- detect_lm(
  cbind(plot_response, 1, plot_x), beta = 1e6,
  cost_adjustment = "BIC", cp_only = FALSE
)
plot_multiple_residual_fit <- plot_residual_fit
plot_multiple_residual_fit@residuals <- cbind(
  plot_residual_fit@residuals, rep(999, length(plot_x))
)
plot_missing_residual_fit <- plot_residual_fit
plot_missing_residual_fit@residuals <- matrix(numeric(), 0, 0)
plot_short_residual_fit <- plot_residual_fit
plot_short_residual_fit@residuals <- plot_residual_fit@residuals[-1, , drop = FALSE]
plot_zero_column_fit <- plot_residual_fit
plot_zero_column_fit@residuals <- matrix(numeric(), length(plot_x), 0)
plot(plot_residual_fit)
plot_regular_residuals <- ggplot2::last_plot()$layers[[3]]$data$y
plot(plot_multiple_residual_fit)
plot_first_residuals <- ggplot2::last_plot()$layers[[3]]$data$y
plot(plot_missing_residual_fit)
plot_missing_residuals <- ggplot2::last_plot()$layers[[3]]$data$y
plot(plot_short_residual_fit)
plot_short_residuals <- ggplot2::last_plot()$layers[[3]]$data$y
plot(plot_zero_column_fit)
plot_zero_column_residuals <- ggplot2::last_plot()$layers[[3]]$data$y
