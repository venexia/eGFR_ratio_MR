make_instrument <- function(
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
  pos,
  p_threshold,
  clump = TRUE,
  clump_kb = 10000,
  clump_r2 = 0.001
) {
  # Load exposure data ----
  message('Load exposure data')

  df <- vroom::vroom(data)

  # Add missing information ----
  message('Add missing information')

  df$Phenotype <- phenotype

  # Format exposure data ----
  message('Format exposure data')

  if (samplesize == "") {
    df <- TwoSampleMR::format_data(
      df,
      type = "exposure",
      phenotype_col = phenotype,
      snp_col = snp,
      effect_allele_col = effect_allele,
      other_allele_col = other_allele,
      eaf_col = eaf,
      beta_col = beta,
      se_col = se,
      pval_col = pval,
      chr_col = chr,
      pos_col = pos
    )
  } else {
    df <- TwoSampleMR::format_data(
      df,
      type = "exposure",
      phenotype_col = phenotype,
      snp_col = snp,
      effect_allele_col = effect_allele,
      other_allele_col = other_allele,
      eaf_col = eaf,
      beta_col = beta,
      se_col = se,
      pval_col = pval,
      samplesize_col = samplesize,
      chr_col = chr,
      pos_col = pos
    )
  }

  # Make instrument ----
  message('Make instrument')

  df <- df[
    df$pval.exposure < p_threshold &
      nchar(df$effect_allele.exposure) == 1 &
      nchar(df$other_allele.exposure) == 1,
  ]

  # Clump instrument ----
  message('Clump instrument')

  if (isTRUE(clump)) {
    df <- TwoSampleMR::clump_data(df, clump_kb = clump_kb, clump_r2 = clump_r2)
  }

  # Return instrument ----
  message('Return instrument')

  return(df)
}
