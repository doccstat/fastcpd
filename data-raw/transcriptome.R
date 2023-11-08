load("data-raw/transcriptome.RData")
transcriptome <- data.frame(ACGH$data, row.names = seq_len(nrow(ACGH$data)))
colnames(transcriptome) <- ACGH$individual
usethis::use_data(transcriptome, overwrite = TRUE)
