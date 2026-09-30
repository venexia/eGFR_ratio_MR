make_instrument_list <- function(
  exposure_files = list.files(
    "data/exposure",
    pattern = ".csv",
    full.names = TRUE
  )
) {
  df <- data.frame(snp = character())

  for (i in exposure_files) {
    tmp <- data.table::fread(i, select = c("SNP"), data.table = FALSE)
    df <- rbind(df, tmp)
  }

  df <- unique(df)
  data.table::fwrite(df, "data/df.csv", row.names = FALSE)

  return(df)
}
