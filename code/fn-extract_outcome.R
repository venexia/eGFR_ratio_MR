extract_outcome <- function(
  instrument,
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
    snps = instrument,
    phenotype_col = phenotype,
    snp_col = snp,
    effect_allele_col = effect_allele,
    other_allele_col = other_allele,
    eaf_col = eaf,
    beta_col = beta(),
    se_col = se,
    pval_col = pval,
    samplesize_col = samplesize,
    chr_col = chr,
    pos_col = pos
  )

  # Return outcome ----
  message('Return outcome')

  return(df)
}
