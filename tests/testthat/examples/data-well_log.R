if (!requireNamespace("ggplot2", quietly = TRUE)) utils::install.packages(
  "ggplot2", repos = "https://cloud.r-project.org", quiet = TRUE
)

result <- fastcpd.mean(well_log, trim = 0.003)
summary(result)
plot(result)
