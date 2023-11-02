if (!requireNamespace("ggplot2", quietly = TRUE)) utils::install.packages(
  "ggplot2", repos = "https://cloud.r-project.org", quiet = TRUE
)

ggplot2::ggplot(bitcoin, ggplot2::aes(x = date, y = price)) +
  ggplot2::geom_line()
