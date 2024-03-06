if (!requireNamespace("ggplot2", quietly = TRUE)) utils::install.packages(
  "ggplot2", repos = "https://cloud.r-project.org", quiet = TRUE
)

result <- fastcpd.mean(
  well_log, beta = (1 + 2) * log(length(well_log)) / 2 * 6, trim = 0.002
)
summary(result)
plot(result)
