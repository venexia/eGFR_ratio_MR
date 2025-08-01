# Load eGFR ratio GWAS ----

df <- vroom::vroom(
  "raw/eGFR11_cys_cre_ratio_imputed.txt.gz",
  col_select = c("SNP", "CHR", "BP", "A1", "ALLELE0")
)

# Duplicate dataste with A1/A0 and A0/A1 ----

df1 <- df
df1 <- dplyr::rename(df1, "allele2" = "A1", "allele1" = "ALLELE0")

df2 <- df
df2 <- dplyr::rename(df2, "allele1" = "A1", "allele2" = "ALLELE0")

df <- rbind(df1, df2)

# Make SNP mappings ----

df$SNP_CAD <- paste0(
  df$CHR,
  ":",
  df$BP,
  "_",
  toupper(df$allele1),
  "_",
  toupper(df$allele2)
)

# Save dataset ----

df <- unique(df[, c("SNP", "SNP_CAD")])
data.table::fwrite(df, "raw/SNPmap.csv")
