# Change-point-only fits retain enough data to plot observations and boundaries.
plot_data <- c(rep(0, 50), rep(4, 50))
plot_cp_only_fit <- detect_mean(
  plot_data, beta = 2, cost_adjustment = "BIC",
  variance_estimation = 1, cp_only = TRUE
)
plot_detail_fit <- detect_mean(
  plot_data, beta = 2, cost_adjustment = "BIC",
  variance_estimation = 1, cp_only = FALSE
)
plot(plot_cp_only_fit)
plot(plot_detail_fit)
