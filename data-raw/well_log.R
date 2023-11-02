well_log <- read.csv("data-raw/well_log.txt", header = FALSE)$V1
well_log <- stats::ts(well_log[well_log > 1e5])
usethis::use_data(well_log, overwrite = TRUE)
