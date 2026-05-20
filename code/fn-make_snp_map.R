make_snp_map <- function() {
  # Load eGFR ratio GWAS ----

  df <- vroom::vroom(
    "raw/eGFR11_cys_cre_ratio_imputed.txt.gz",
    col_select = c("SNP", "CHR", "BP", "A1", "ALLELE0")
  )

  df <- dplyr::rename(df, "chr" = "CHR", "pos" = "BP")

  # Duplicate dataste with A1/A0 and A0/A1 ----

  df1 <- df
  df1 <- dplyr::rename(df1, "effect_allele" = "A1", "other_allele" = "ALLELE0")

  df2 <- df
  df2 <- dplyr::rename(df2, "other_allele" = "A1", "effect_allele" = "ALLELE0")

  df <- rbind(df1, df2)

  # Filter to biallelic SNPs ----

  df <- df[nchar(df$other_allele) == 1 & nchar(df$effect_allele) == 1, ]

  # Save dataset ----

  data.table::fwrite(df, "raw/SNPmap.csv")

  # Return dataset ----

  return(SNPmap)
}
