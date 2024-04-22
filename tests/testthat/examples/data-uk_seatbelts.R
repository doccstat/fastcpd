if (
  requireNamespace("ggplot2", quietly = TRUE) &&
    requireNamespace("lubridate", quietly = TRUE) &&
    requireNamespace("zoo", quietly = TRUE)
) {
  result_ar <- fastcpd.ar(diff(uk_seatbelts[, "drivers"], 12), 1, beta = "BIC")
  summary(result_ar)
  plot(result_ar)

  result_lm <- suppressMessages(fastcpd.lm(
    diff(uk_seatbelts[, c("drivers", "kms", "PetrolPrice", "law")], lag = 12)
  ))
  cp_dates <- as.Date("1969-01-01", format = "%Y-%m-%d")
  cp_dates <- cp_dates + lubridate::period(month = 1 + result_lm@cp_set + 12)
  cp_dates <- zoo::as.yearmon(cp_dates)

  dates <- zoo::as.yearmon(time(uk_seatbelts))

  uk_seatbelts_df <- data.frame(
    dates = dates,
    drivers = c(uk_seatbelts[, "drivers"]),
    color = as.factor((dates < cp_dates[1]) + (dates < cp_dates[2]))
  )

  ggplot2::ggplot() +
    ggplot2::geom_line(
      data = uk_seatbelts_df,
      mapping = ggplot2::aes(x = dates, y = drivers, color = color)
    ) +
    ggplot2::geom_vline(
      xintercept = cp_dates,
      linetype = "dashed",
      color = "red"
    ) +
    zoo::scale_x_yearmon() +
    ggplot2::annotate(
      "text",
      x = cp_dates,
      y = 1025,
      label = as.character(cp_dates),
      color = "blue"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none")
}
