hpse_path <- "raw/HPSE_Q9Y251_OID30409_v1_Cardiometabolic_II/"
hpse_files <- list.files(
  hpse_path,
  pattern = "HPSE:Q9Y251:OID30409:v1:Cardiometabolic_II.gz",
  full.names = TRUE
)

df <- data.table::fread(hpse_files[1], data.table = FALSE, nrows = 0)

for (i in hpse_files) {
  tmp <- data.table::fread(i, data.table = FALSE)
  df <- rbind(df, tmp)
}

readr::write_rds(df, "data/hpse.rds", compress = "gz")
