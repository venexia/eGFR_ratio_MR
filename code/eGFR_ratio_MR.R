# Source functions ----

lapply(
  list.files("code", full.names = TRUE, pattern = "fn-"),
  source
)

# Make instrument ----

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

out <- extract_outcome(
  data = "",
  phenotype = "",
  samplesize = NA,
  snps = exp$SNP
)

# Harmonize data ----

dat <- TwoSampleMR::harmonise_data(exposure_dat = exp, outcome_dat = out)

# Perform MR ----

mr <- TwoSampleMR::mr(dat)

# Record extra information ----

mr$nsnp.exposure <- nrow(exp)
