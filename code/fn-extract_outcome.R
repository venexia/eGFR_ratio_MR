extract_outcome <- function(data, phenotype, samplesize, snps) {
  # Load outcome data ----
  df <- vroom::vroom(data)

  # Add missing information ----

  df$Phenotype <- phenotype
  df$samplesize <- as.numeric(samplesize)

  df <- TwoSampleMR::format_data(
    out,
    type = "outcome",
    snps = snps,
    phenotype_col = "Phenotype",
    snp_col = "RSID",
    effect_allele_col = "Allele1",
    other_allele_col = "Allele0",
    eaf_col = "Freq1",
    beta_col = "Effect",
    se_col = "StdErr",
    pval_col = "P-value",
    samplesize_col = "n_total_sum",
    chr_col = "Chr",
    pos_col = "Pos_b37"
  )

  # Return outcome ----

  return(df)
}
