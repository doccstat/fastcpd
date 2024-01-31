for (package in c("ggplot2", "lubridate", "zoo")) {
  if (!requireNamespace(package, quietly = TRUE)) utils::install.packages(
    package, repos = "https://cloud.r-project.org", quiet = TRUE
  )
}

result_ar <- fastcpd.ar(diff(uk_seatbelts[, "drivers"], 12), 1, beta = "BIC")
summary(result_ar)
plot(result_ar)

result_lm <- suppressMessages(fastcpd.lm(
  diff(uk_seatbelts[, c("drivers", "kms", "PetrolPrice", "law")], lag = 12)
))
cp_dates <- as.Date("1969-01-01", format = "%Y-%m-%d")
cp_dates <- cp_dates + lubridate::period(month = 1 + result_lm@cp_set + 12)
cp_dates <- zoo::as.yearmon(cp_dates)

uk_seatbelts_df <- data.frame(
  dates = zoo::as.yearmon(time(uk_seatbelts)),
  drivers = c(uk_seatbelts[, "drivers"])
)

ggplot2::ggplot() +
  ggplot2::geom_line(
    data = uk_seatbelts_df,
    mapping = ggplot2::aes(x = dates, y = drivers)
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
  )
