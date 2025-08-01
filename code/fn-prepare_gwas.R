library(rlang)
library(magrittr)

prepare_gwas <- function(
  map_snps,
  type,
  data,
  phenotype,
  phenotype_short,
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
  id,
  p_threshold, # only needed when type = "exposure"
  clump, # only needed when type = "exposure"
  clump_kb, # only needed when type = "exposure"
  clump_r2, # only needed when type = "exposure"
  instrument # only needed when type = "outcome"
) {
  if (id != "" & type == "exposure") {
    df <- TwoSampleMR::extract_instruments(outcomes = id)
  } else if (id != "" & type == "outcome") {
    df <- batch_extract_outcome(snp_list = instrument, outcomes = id)
  } else {
    # Load GWAS data ----
    message('Load GWAS data')

    df <- vroom::vroom(data)

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

    # Map SNPs ----

    if (map_snps == "TRUE") {
      # Load SNP map ----
      message('Load SNP map')

      map <- vroom::vroom(
        "raw/SNPmap.csv",
        col_select = c("SNP", paste0("SNP_", phenotype_short))
      )

      message(paste0("Map contains "), nrow(map), " SNPs")

      # Restrict SNP map to SNPs in GWAS ----
      message('Restrict SNP map to SNPs in GWAS')

      df$SNP <- toupper(df$SNP)
      map <- map[map$SNP %in% df$SNP, ]

      message(paste0("Map contains "), nrow(map), " SNPs")

      if (type == "outcome") {
        # Restrict SNP map to SNPs in instrument ----
        message('Restrict SNP map to SNPs in instrument')
        map <- map[map$SNP %in% instrument, ]
        message(paste0("Map contains "), nrow(map), " SNPs")
      }

      # Label SNP map ----
      message('Label SNP map')

      df <- df %>%
        dplyr::rename(!!paste0("SNP_", phenotype_short) := SNP)

      # Map SNPs ----
      message('Map SNPs')
      df <- merge(df, map, by = paste0("SNP_", phenotype_short), all.x = TRUE)
    }

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
  }

  # Standardize betas and SEs for lifetime smoking index

  if (phenotype == "Lifetime smoking index") {
    message('Standardize betas and SEs for lifetime smoking index')
    LSI_SD <- 0.6940093

    if (type == "exposure") {
      df$beta.exposure <- df$beta.exposure / LSI_SD
      df$se.exposure <- df$se.exposure / LSI_SD
    }

    if (type == "outcome") {
      df$beta.outcome <- df$beta.outcome / LSI_SD
      df$se.outcome <- df$se.outcome / LSI_SD
    }
  }

  # Label phenotype ----
  message('Label phenotype')

  if (type == "exposure") {
    df$exposure <- phenotype
  }

  if (type == "outcome") {
    df$outcome <- phenotype
  }

  # Return gwas ----
  message('Return gwas')

  return(df)
}
