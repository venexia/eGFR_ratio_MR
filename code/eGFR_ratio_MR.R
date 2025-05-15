# Source functions ----
message("Source functions")

lapply(
  list.files("code", full.names = TRUE, pattern = "fn-"),
  source
)

# Make GWAS list ----
message("Make GWAS list")

source("code/gwas.R")

# Make empty results dataframe ----
message("Make empty results dataframe")

results <- NULL

# Prepare exposure GWAS ----
message("Prepare exposure GWAS")

exp <- NULL
exposures <- gwas[gwas$type == "exposure", ]
exposures <- exposures[1, ]

for (i in 1:nrow(exposures)) {
  message(paste0("Preparing ", exposures[i, "phenotype"]))

  tmp <- prepare_gwas(
    type = "exposure",
    data = exposures[i, "data"],
    phenotype = exposures[i, "phenotype"],
    snp_col = exposures[i, "snp_col"],
    effect_allele_col = exposures[i, "effect_allele_col"],
    other_allele_col = exposures[i, "other_allele_col"],
    eaf_col = exposures[i, "eaf_col"],
    beta_col = exposures[i, "beta_col"],
    se_col = exposures[i, "se_col"],
    pval_col = exposures[i, "pval_col"],
    samplesize_col = exposures[i, "samplesize_col"],
    samplesize = exposures[i, "samplesize"],
    chr_col = exposures[i, "chr_col"],
    pos_col = exposures[i, "pos_col"],
    p_threshold = 5e-8,
    clump = TRUE,
    clump_kb = 10000,
    clump_r2 = 0.001
  )

  exp <- rbind(exp, tmp)
}

# Prepare outcome GWAS ----
message("Prepare outcome GWAS")

out <- NULL
outcomes <- gwas[gwas$type == "outcome", ]

for (i in 1:nrow(outcomes)) {
  message(paste0("Preparing "), outcomes[i, "phenotype"])

  tmp <- prepare_gwas(
    type = "outcome",
    data = outcomes[i, "data"],
    phenotype = outcomes[i, "phenotype"],
    snp_col = outcomes[i, "snp_col"],
    effect_allele_col = outcomes[i, "effect_allele_col"],
    other_allele_col = outcomes[i, "other_allele_col"],
    eaf_col = outcomes[i, "eaf_col"],
    beta_col = outcomes[i, "beta_col"],
    se_col = outcomes[i, "se_col"],
    pval_col = outcomes[i, "pval_col"],
    samplesize_col = outcomes[i, "samplesize_col"],
    samplesize = outcomes[i, "samplesize"],
    chr_col = outcomes[i, "chr_col"],
    pos_col = outcomes[i, "pos_col"],
    instrument = exp$SNP
  )

  out <- rbind(out, tmp)
}

# Harmonize data ----
message("Harmonize data")

dat <- TwoSampleMR::harmonise_data(exposure_dat = exp, outcome_dat = out)

# Perform MR ----
message("Perform MR")

mr <- TwoSampleMR::mr(dat)

# Record extra information ----
message("Record extra information")

mr$nsnp.exposure <- nrow(exp)

# Add to results dataframe ----
message("Add to results dataframe")

results <- rbind(results, mr)

# Format results dataframe ----
message("Format results dataframe")

results$lci <- exp(results$b - qnorm(0.975) * results$se)
results$uci <- exp(results$b + qnorm(0.975) * results$se)
results$or <- exp(results$b)

results <- results[, c(
  "exposure",
  "outcome",
  "method",
  "nsnp.exposure",
  "nsnp",
  "or",
  "lci",
  "uci",
  "pval"
)]

# Save results dataframe ----
message("Save results dataframe")

data.table::fwrite(results, "output/results.csv", row.names = FALSE)
