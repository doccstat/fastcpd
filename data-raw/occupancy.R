occupancy <- read.csv("data-raw/occupancy.txt", header = TRUE)
usethis::use_data(occupancy, overwrite = TRUE)
