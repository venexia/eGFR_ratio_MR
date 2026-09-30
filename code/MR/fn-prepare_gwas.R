library(rlang)
library(magrittr)

prepare_gwas <- function(
  gwas,
  type,
  p_threshold, # only needed when type = "exposure"
  clump, # only needed when type = "exposure"
  clump_kb, # only needed when type = "exposure"
  clump_r2, # only needed when type = "exposure"
  instrument, # only needed when type = "outcome"
  save = TRUE
) {
  # Check if file already exists ----
  gwas_filename <- paste0("data/", type, "/", gwas$phenotype_short, ".csv")

  if (file.exists(gwas_filename)) {
    message("Skipped. File already exists.")
    df <- data.table::fread(gwas_filename, data.table = FALSE)
  } else {
    if (gwas$id != "") {
      # Extract from IEU Open GWAS if present ----
      if (type == "exposure") {
        df <- TwoSampleMR::extract_instruments(outcomes = gwas$id)
      } else if (type == "outcome") {
        df <- batch_extract_outcome(snp_list = instrument, outcomes = gwas$id)
      }
      df[, colnames(df)[
        grepl("samplesize", colnames(df)) |
          grepl("ncase", colnames(df)) |
          grepl("ncontrol", colnames(df))
      ]] <- NULL
    } else {
      # Load GWAS data ----
      df <- vroom::vroom(gwas$data, show_col_types = FALSE)

      # Format GWAS data ----
      df <- dplyr::rename(
        df,
        effect_allele = all_of(gwas$effect_allele_col),
        other_allele = all_of(gwas$other_allele_col),
        eaf = all_of(gwas$eaf_col),
        beta = all_of(gwas$beta_col),
        se = all_of(gwas$se_col),
        pval = all_of(gwas$pval_col),
        chr = all_of(gwas$chr_col),
        pos = all_of(gwas$pos_col)
      )

      # Format SNP col -----
      if (gwas$snp_col != "") {
        df <- dplyr::rename(
          df,
          SNP = all_of(gwas$snp_col),
        )
      }

      # Convert to numeric ----
      cols <- c("eaf", "beta", "se", "pval")
      df[cols] <- lapply(df[cols], as.numeric)

      # If applicable: map SNPs ----

      if (gwas$map_snps != "") {
        ## Identify map variables ----
        map_snps <- strsplit(gwas$map_snps, ";")[[1]]

        ## Load SNP map ----
        map <- vroom::vroom(
          "raw/SNPmap.csv"
        )

        if (type == "outcome") {
          ## Restrict SNP map to SNPs in instrument ----
          message('Restrict SNP map to SNPs in instrument')
          map <- map[map$SNP %in% instrument, ]
          message(paste0("Map contains "), nrow(map), " SNPs")
        }

        n_orig <- nrow(df)

        ## Map SNPs ----
        df <- merge(df, map, by = map_snps)
        message(paste0(nrow(df), " / ", n_orig, " mapped!"))
      }

      # Reformat data ----
      df <- TwoSampleMR::format_data(
        df,
        type = type,
        snp_col = "SNP",
        effect_allele_col = "effect_allele",
        other_allele_col = "other_allele",
        eaf_col = "eaf",
        beta_col = "beta",
        se_col = "se",
        pval_col = "pval",
        chr_col = "chr",
        pos_col = "pos"
      )

      # If applicable: exposure specific preparations ----
      if (type == "exposure") {
        # Restrict to bialleic SNPs with known effect estimates and p-values ----
        df <- df[
          !is.na(df$pval.exposure) &
            !is.na(df$beta.exposure) &
            nchar(df$effect_allele.exposure) == 1 &
            nchar(df$other_allele.exposure) == 1,
        ]

        if (sum(df$pval.exposure < 5e-8) > 0) {
          # Restrict to genome-wide significant SNPs ----

          df <- df[
            df$pval.exposure < p_threshold,
          ]

          if (nrow(df) > 0) {
            # Clump instrument ----
            if (isTRUE(clump)) {
              df <- TwoSampleMR::clump_data(
                df,
                clump_kb = clump_kb,
                clump_r2 = clump_r2
              )
            }
          }
        } else {
          message("No genome-wide significant SNPs")
          df <- data.frame()
        }
      }

      # If applicable: outcome specific preparations ----
      if (type == "outcome") {
        # Filter to instrument SNPs ----
        df <- df[df$SNP %in% instrument, ]
      }
    }

    if (nrow(df) > 0) {
      # If applicable: standardize betas and SEs ---- ----

      SD <- as.numeric(gwas$SD)

      if (!is.na(SD)) {
        message('Standardize betas and SEs')

        if (type == "exposure") {
          df$beta.exposure <- df$beta.exposure / SD
          df$se.exposure <- df$se.exposure / SD
        }

        if (type == "outcome") {
          df$beta.outcome <- df$beta.outcome / SD
          df$se.outcome <- df$se.outcome / SD
        }
      }

      # Label phenotype ---- ----
      if (type == "exposure" & nrow(df) > 0) {
        df$exposure <- gwas$phenotype
      }

      if (type == "outcome" & nrow(df) > 0) {
        df$outcome <- gwas$phenotype
      }

      # Add samplesizes ---- ----

      samplesize_col <- paste0("samplesize.", type)
      if (!(samplesize_col %in% colnames(df))) {
        message("Add ", samplesize_col)
        df[[samplesize_col]] <- gwas$samplesize
      }

      if (gwas$category == "binary") {
        ncase_col <- paste0("ncase.", type)
        ncontrol_col <- paste0("ncontrol.", type)

        if (!(ncase_col %in% colnames(df))) {
          message("Add ", ncase_col)
          df[[ncase_col]] <- gwas$ncase
        }

        if (!(ncontrol_col %in% colnames(df))) {
          message("Add ", ncontrol_col)
          df[[ncontrol_col]] <- gwas$ncontrol
        }
      }
    }
  }

  # Save gwas ----

  if (!is.null(nrow(df)) & isTRUE(save)) {
    message('Save gwas')
    data.table::fwrite(
      df,
      gwas_filename,
      row.names = FALSE
    )
  }

  # Return gwas ----
  message('Return gwas')

  return(df)
}
