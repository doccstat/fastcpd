well_log <- read.csv("data-raw/well_log.txt", header = FALSE)$V1
well_log <- stats::ts(well_log)
usethis::use_data(well_log, overwrite = TRUE)
