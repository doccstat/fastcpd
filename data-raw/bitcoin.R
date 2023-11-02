bitcoin_json <- rjson::fromJSON(file = "data-raw/bitcoin.json")
bitcoin_ <- data.frame(do.call(Map, c(f = rbind, bitcoin_json$`market-price`)))
bitcoin <- data.frame(
  date = as.POSIXct(bitcoin_$x / 1000, origin = "1969-12-31", tz = "UTC"),
  price = bitcoin_$y
)
usethis::use_data(bitcoin, overwrite = TRUE)
