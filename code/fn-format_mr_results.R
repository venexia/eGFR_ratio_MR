format_mr_results <- function(filepath = "data/mr/", metadata) {
  
  # Get files ----
  message("Get files")
  
  mr_files <- list.files(
    path = filepath,
    pattern = "\\.csv$",
    full.names = TRUE
  )
  
  results <- data.table::rbindlist(
    lapply(mr_files, data.table::fread),
    use.names = TRUE,
    fill = TRUE
  )
  
  
  # Add GWAS meta data ----
  message("Add GWAS meta data")
  
  results <- merge(
    results,
    metadata[, c("phenotype", "category", "ukb")],
    by.x = "outcome",
    by.y = "phenotype",
    all.x = TRUE
  )
  
  results <- dplyr::rename(results, "ukb.outcome" = "ukb")
  
  results <- merge(
    results,
    metadata[, c("phenotype", "ukb")],
    by.x = "exposure",
    by.y = "phenotype",
    all.x = TRUE
  )
  
  results <- dplyr::rename(results, "ukb.exposure" = "ukb")
  
  # Format results dataframe ----
  message("Format results dataframe")
  
  results$b_lci <- results$b - qnorm(0.975) * results$se
  results$b_uci <- results$b + qnorm(0.975) * results$se
  results$or_lci <- exp(results$b_lci)
  results$or_uci <- exp(results$b_uci)
  results$or <- exp(results$b)
  
  results$est <- ifelse(results$category == "continuous", results$b, results$or)
  results$lci <- ifelse(
    results$category == "continuous",
    results$b_lci,
    results$or_lci
  )
  results$uci <- ifelse(
    results$category == "continuous",
    results$b_uci,
    results$or_uci
  )
  results$category <- ifelse(
    results$category == "continuous",
    results$category,
    "binary"
  )
  
  results <- results[, c(
    "exposure",
    "outcome",
    "method",
    "nsnp.exposure",
    "nsnp",
    "b",
    "b_lci",
    "b_uci",
    "category",
    "est",
    "lci",
    "uci",
    "pval"
  )]
  
}