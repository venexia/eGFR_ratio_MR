extract_outcome <- function(
  snps,
  data,
  phenotype,
  snp,
  effect_allele,
  other_allele,
  eaf,
  beta,
  se,
  pval,
  samplesize,
  chr,
  pos
) {
  # Load outcome data ----
  message('Load outcome data')
  df <- vroom::vroom(data)

  # Add missing information ----
  message('Add missing information')

  df$Phenotype <- phenotype

  # Format data as an outcome ----
  message('Format data as an outcome')

  df <- TwoSampleMR::format_data(
    df,
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
  message('Return outcome')

  return(df)
}
