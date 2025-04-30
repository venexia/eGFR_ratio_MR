# Source functions ----
message("Source functions")

lapply(
  list.files("code", full.names = TRUE, pattern = "fn-"),
  source
)

# Make outcomes list ----
message("Make outcomes list")

source("code/outcomes.R")

# Run make instrument function ----
message("Run make instrument function")

exp <- make_instrument(
  data = "raw/eGFR11_cys_cre_ratio_imputed.txt.gz",
  phenotype = "eGFRcys_eGFRcr",
  samplesize = NA,
  p_threshold = 5e-8,
  clump = TRUE,
  clump_kb = 10000,
  clump_r2 = 0.001
)

# Extract outcome data ----
message("Extract outcome data")

for (i in 1:nrow(outcomes))
  out <- extract_outcome(
    snps = exp$SNP,
    data = outcomes[i, "path"],
    phenotype = outcomes[i, "phenotype"],
    snp = outcomes[i, "snp"],
    effect_allele = outcomes[i, "effect_allele"],
    other_allele = outcomes[i, "other_allele"],
    eaf = outcomes[i, "eaf"],
    beta = outcomes[i, "beta"],
    se = outcomes[i, "se"],
    pval = outcomes[i, "pval"],
    samplesize = outcomes[i, "samplesize"],
    chr = outcomes[i, "chr"],
    pos = outcomes[i, "pos"]
  )

# Harmonize data ----
message("Harmonize data")

dat <- TwoSampleMR::harmonise_data(exposure_dat = exp, outcome_dat = out)

# Perform MR ----
message("Perform MR")

mr <- TwoSampleMR::mr(dat)

# Record extra information ----
message("Record extra information")

mr$nsnp.exposure <- nrow(exp)
