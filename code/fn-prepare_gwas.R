prepare_gwas <- function(
  type,
  data,
  phenotype,
  snp_col,
  effect_allele_col,
  other_allele_col,
  eaf_col,
  beta_col,
  se_col,
  pval_col,
  samplesize_col,
  samplesize,
  chr_col,
  pos_col,
  p_threshold, # only needed when type = "exposure"
  clump, # only needed when type = "exposure"
  clump_kb, # only needed when type = "exposure"
  clump_r2, # only needed when type = "exposure"
  instrument # only needed when type = "outcome"
) {
  # Load GWAS data ----
  message('Load GWAS data')

  df <- vroom::vroom(data)

  # Label phenotype ----
  message('Label phenotype')

  df$phenotype <- phenotype

  # Add sample size ----

  if (samplesize_col == "") {
    message('Add sample size')
    df$samplesize <- samplesize
    samplesize_col <- "samplesize"
  }

  # Format data ----
  message('Format data')

  df <- TwoSampleMR::format_data(
    df,
    type = type,
    phenotype_col = "phenotype",
    snp_col = snp_col,
    effect_allele_col = effect_allele_col,
    other_allele_col = other_allele_col,
    eaf_col = eaf_col,
    beta_col = beta_col,
    se_col = se_col,
    pval_col = pval_col,
    samplesize_col = samplesize_col,
    chr_col = chr_col,
    pos_col = pos_col
  )

  # Exposure specific preparations ----

  if (type == "exposure") {
    # Make instrument ----
    message('Make instrument')

    df <- df[
      df$pval.exposure < p_threshold &
        !is.na(df$pval.exposure) &
        !is.na(df$beta.exposure) &
        nchar(df$effect_allele.exposure) == 1 &
        nchar(df$other_allele.exposure) == 1,
    ]

    # Clump instrument ----
    message('Clump instrument')

    if (isTRUE(clump)) {
      df <- TwoSampleMR::clump_data(
        df,
        clump_kb = clump_kb,
        clump_r2 = clump_r2
      )
    }
  }

  # Outcome specific preparations ----

  if (type == "outcome") {
    # Filter to instrument SNPs ----
    message('Filter to instrument SNPs')
    df <- df[df$SNP %in% instrument, ]
  }

  # Return gwas ----
  message('Return gwas')

  return(df)
}
