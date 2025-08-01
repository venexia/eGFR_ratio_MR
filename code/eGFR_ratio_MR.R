# Make GWAS list ----
message("Make GWAS list")

source("code/gwas.R")
rm(list = setdiff(ls(), "gwas"))

# Source functions ----
message("Source functions")

lapply(
  list.files("code", full.names = TRUE, pattern = "fn-"),
  source
)

# Make SNP map ----
message("Make SNP map")

# source("code/make_SNPmap.R")

# Make empty results dataframe ----
message("Make empty results dataframe")

results <- NULL

# # Prepare exposure GWAS ----
# message("Prepare exposure GWAS")
#
# exp <- data.frame(
#   SNP = character(),
#   chr.exposure = character(),
#   pos.exposure = character(),
#   effect_allele.exposure = character(),
#   other_allele.exposure = character(),
#   eaf.exposure = character(),
#   beta.exposure = character(),
#   se.exposure = character(),
#   pval.exposure = character(),
#   samplesize.exposure = character(),
#   exposure = character(),
#   mr_keep.exposure = character(),
#   pval_origin.exposure = character()
# )
#
# for (i in 1:nrow(gwas)) {
#   message(paste0(
#     "Preparing exposure ",
#     i,
#     " of ",
#     nrow(gwas),
#     ": ",
#     gwas[i, "phenotype"]
#   ))
#
#   tmp <- prepare_gwas(
#     map_snps = gwas[i, "map_snps"],
#     type = "exposure",
#     data = gwas[i, "data"],
#     phenotype = gwas[i, "phenotype"],
#     phenotype_short = gwas[i, "phenotype_short"],
#     snp_col = gwas[i, "snp_col"],
#     effect_allele_col = gwas[i, "effect_allele_col"],
#     other_allele_col = gwas[i, "other_allele_col"],
#     eaf_col = gwas[i, "eaf_col"],
#     beta_col = gwas[i, "beta_col"],
#     se_col = gwas[i, "se_col"],
#     pval_col = gwas[i, "pval_col"],
#     samplesize_col = gwas[i, "samplesize_col"],
#     samplesize = gwas[i, "samplesize"],
#     chr_col = gwas[i, "chr_col"],
#     pos_col = gwas[i, "pos_col"],
#     id = gwas[i, "id"],
#     p_threshold = 5e-8,
#     clump = TRUE,
#     clump_kb = 10000,
#     clump_r2 = 0.001
#   )
#
#   if (!is.null(nrow(tmp))) {
#     exp <- plyr::rbind.fill(exp, tmp)
#   }
# }
#
# data.table::fwrite(exp, "data/exposures.csv", row.names = FALSE)

exp <- data.table::fread("data/exposures.csv")

# Prepare outcome GWAS ----
message("Prepare outcome GWAS")

out <- NULL

for (i in 1:nrow(gwas)) {
  message(paste0(
    "Preparing outcome ",
    i,
    " of ",
    nrow(gwas),
    ": ",
    gwas[i, "phenotype"]
  ))

  tmp <- prepare_gwas(
    map_snps = gwas[i, "map_snps"],
    type = "outcome",
    data = gwas[i, "data"],
    phenotype = gwas[i, "phenotype"],
    phenotype_short = gwas[i, "phenotype_short"],
    snp_col = gwas[i, "snp_col"],
    effect_allele_col = gwas[i, "effect_allele_col"],
    other_allele_col = gwas[i, "other_allele_col"],
    eaf_col = gwas[i, "eaf_col"],
    beta_col = gwas[i, "beta_col"],
    se_col = gwas[i, "se_col"],
    pval_col = gwas[i, "pval_col"],
    samplesize_col = gwas[i, "samplesize_col"],
    samplesize = gwas[i, "samplesize"],
    chr_col = gwas[i, "chr_col"],
    pos_col = gwas[i, "pos_col"],
    id = gwas[i, "id"],
    instrument = unique(exp$SNP)
  )

  if (!is.null(nrow(tmp))) {
    out <- plyr::rbind.fill(out, tmp)
  }
}

# Harmonize data ----
message("Harmonize data")

dat <- TwoSampleMR::harmonise_data(exposure_dat = exp, outcome_dat = out)

# Perform MR ----
message("Perform MR")

mr <- TwoSampleMR::mr(dat)

# Record extra information ----
message("Record extra information")

tmp <- as.data.frame(table(exp$exposure))
tmp <- dplyr::rename(tmp, "exposure" = "Var1", "nsnp.exposure" = "Freq")
tmp$exposure <- as.character(tmp$exposure)
mr <- merge(mr, tmp, by = "exposure", all.x = TRUE)

# Add to results dataframe ----
message("Add to results dataframe")

results <- rbind(results, mr)

# Add GWAS meta data ----
message("Add GWAS meta data")

results <- merge(
  results,
  gwas[, c("phenotype", "category", "ukb")],
  by.x = "outcome",
  by.y = "phenotype",
  all.x = TRUE
)

# Format results dataframe ----
message("Format results dataframe")

results$b_lci <- results$b - qnorm(0.975) * results$se
results$b_uci <- results$b + qnorm(0.975) * results$se
results$or_lci <- exp(results$lci)
results$or_uci <- exp(results$uci)
results$or <- exp(results$b)

results$est <- ifelse(results)

results <- results[, c(
  "exposure",
  "outcome",
  "method",
  "nsnp.exposure",
  "nsnp",
  "b",
  "lci",
  "uci",
  "pval"
)]

# Save results dataframe ----
message("Save results dataframe")

data.table::fwrite(results, "output/results.csv", row.names = FALSE)
