library(rlang)
library(magrittr)

prepare_gwas <- function(
  gwas,
  type,
  p_threshold, # only needed when type = "exposure"
  clump, # only needed when type = "exposure"
  clump_kb, # only needed when type = "exposure"
  clump_r2, # only needed when type = "exposure"
  instrument # only needed when type = "outcome"
) {
  # Check if file already exists
  message('Check if file already exists')
  gwas_filename <- paste0("data/", type, "/", gwas$phenotype_short, ".csv")

  if (file.exists(gwas_filename)) {
    message("Skipped. File already exists.")
    df <- data.table::fread(gwas_filename, data.table = FALSE)
  } else {
    # Determine whether data is in IEU Open GWAS
    message('Determine whether data is in IEU Open GWAS')

    if (gwas$id != "") {
      # Extract from IEU Open GWAS
      message('Extract from IEU Open GWAS')
      if (type == "exposure") {
        df <- TwoSampleMR::extract_instruments(outcomes = gwas$id)
      } else if (type == "outcome") {
        df <- batch_extract_outcome(snp_list = instrument, outcomes = gwas$id)
      }
    } else {
      # Load GWAS data ----
      message('Load GWAS data')

      df <- vroom::vroom(gwas$data, show_col_types = FALSE)

      # Add sample sizes ----

      if (gwas$samplesize_col == "") {
        message('Add sample size')
        df$samplesize <- gwas$samplesize
        samplesize_col <- "samplesize"
      }

      if (gwas$category == "binary") {
        message('Add ncase')
        df$ncase <- gwas$ncase
        ncase_col <- "ncase"
      }

      if (gwas$category == "binary") {
        message('Add ncontrol')
        df$ncontrol <- gwas$ncontrol
        ncontrol_col <- "ncontrol"
      }

      # Format GWAS data ----
      message('Format GWAS data')

      df <- dplyr::rename(
        df,
        effect_allele = all_of(gwas$effect_allele_col),
        other_allele = all_of(gwas$other_allele_col),
        eaf = all_of(gwas$eaf_col),
        beta = all_of(gwas$beta_col),
        se = all_of(gwas$se_col),
        pval = all_of(gwas$pval_col),
        samplesize = all_of(gwas$samplesize_col),
        chr = all_of(gwas$chr_col),
        pos = all_of(gwas$pos_col)
      )

      if (gwas$category == "binary") {
        df <- dplyr::rename(
          df,
          ncase = all_of(ncase_col),
          ncontrol = all_of(ncontrol_col)
        )
      }

      if (gwas$snp_col != "") {
        df <- dplyr::rename(
          df,
          SNP = all_of(gwas$snp_col),
        )
      }

      # Convert to numeric
      cols <- c("eaf", "beta", "se", "pval", "samplesize")
      if (gwas$category == "binary") {
        cols <- c(cols, "ncase", "ncontrol")
      }
      df[cols] <- lapply(df[cols], as.numeric)

      # Map SNPs ----

      if (gwas$map_snps != "") {
        # Identify map variables ----
        message('Identify map variables')

        map_snps <- strsplit(gwas$map_snps, ";")[[1]]

        # Load SNP map ----
        message('Load SNP map')

        map <- vroom::vroom(
          "raw/SNPmap.csv"
        )

        if (type == "outcome") {
          # Restrict SNP map to SNPs in instrument ----
          message('Restrict SNP map to SNPs in instrument')
          map <- map[map$SNP %in% instrument, ]
          message(paste0("Map contains "), nrow(map), " SNPs")
        }

        n_orig <- nrow(df)

        # Map SNPs ----
        message('Map SNPs')
        df <- merge(df, map, by = map_snps)
        message(paste0(nrow(df), " / ", n_orig, " mapped!"))
      }

      # Reformat data ----
      message('Reformat data')

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
        samplesize_col = "samplesize",
        ncase = "ncase",
        ncontrol = "ncontrol",
        chr_col = "chr",
        pos_col = "pos"
      )

      # Exposure specific preparations ----

      if (type == "exposure") {
        # Check necessary instrument info is present ----
        message('Check necessary instrument info is present')

        df <- df[
          df$pval.exposure < p_threshold &
            !is.na(df$pval.exposure) &
            !is.na(df$beta.exposure) &
            nchar(df$effect_allele.exposure) == 1 &
            nchar(df$other_allele.exposure) == 1,
        ]

        if (nrow(df) > 0) {
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
      }

      # Outcome specific preparations ----

      if (type == "outcome") {
        # Filter to instrument SNPs ----
        message('Filter to instrument SNPs')
        df <- df[df$SNP %in% instrument, ]
      }
    }

    # Standardize betas and SEs

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

    # Label phenotype ----
    message('Label phenotype')

    if (type == "exposure" & nrow(df) > 0) {
      df$exposure <- gwas$phenotype
    }

    if (type == "outcome" & nrow(df) > 0) {
      df$outcome <- gwas$phenotype
    }

    # Remove irrelevant columns ----
    message('Remove irrelevant columns')

    if (gwas$category == "continuous") {
      rm_cols <- intersect(
        colnames(df),
        c(
          "ncase.exposure",
          "ncontrol.exposure",
          "ncase.outcome",
          "ncontrol.outcome"
        )
      )
      df[, rm_cols] <- NULL
    }

    # Save gwas ----
    message('Save gwas')

    if (!is.null(nrow(df))) {
      data.table::fwrite(
        df,
        gwas_filename,
        row.names = FALSE
      )
    }
  }

  # Return gwas ----
  message('Return gwas')

  return(df)
}
