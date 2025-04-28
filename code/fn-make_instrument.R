make_instrument <- function(
  data,
  phenotype,
  samplesize,
  p_threshold,
  clump = TRUE,
  clump_kb = 10000,
  clump_r2 = 0.001
) {
  # Load exposure data ----

  df <- vroom::vroom(data)

  # Add missing information ----

  df$Phenotype <- phenotype
  df$samplesize <- as.numeric(samplesize)

  # Format exposure data ----

  df <- TwoSampleMR::format_data(
    df,
    type = "exposure",
    phenotype_col = "Phenotype",
    snp_col = "SNP",
    effect_allele_col = "A1",
    other_allele_col = "ALLELE0",
    eaf_col = "A1FREQ",
    beta_col = "BETA",
    se_col = "SE",
    pval_col = "P_BOLT_LMM_INF",
    samplesize_col = "samplesize",
    chr_col = "CHR",
    pos_col = "BP"
  )

  # Make instrument ----

  df <- df[
    df$pval.exposure <= p_threshold &
      nchar(df$effect_allele.exposure) == 1 &
      nchar(df$other_allele.exposure) == 1,
  ]

  # Clump instrument ----

  if (isTRUE(clump)) {
    df <- TwoSampleMR::clump_data(df, clump_kb = 10000, clump_r2 = 0.001)
  }

  # Return instrument ----

  return(df)
}
