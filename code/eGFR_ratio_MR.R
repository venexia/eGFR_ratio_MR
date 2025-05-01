# Source functions ----
message("Source functions")

lapply(
  list.files("code", full.names = TRUE, pattern = "fn-"),
  source
)

# Make GWAS list ----
message("Make GWAS list")

source("code/gwas.R")

# Run make instrument function ----
message("Run make instrument function")

exposures <- gwas[gwas$type == "exposure" & gwas$phenotype == "eGFR_ratio", ]

for (i in 1:nrow(exposures)) {
  exp <- make_instrument(
    data = exposures[i, "data"],
    phenotype = exposures[i, "phenotype"],
    snp = exposures[i, "snp"],
    effect_allele = exposures[i, "effect_allele"],
    other_allele = exposures[i, "other_allele"],
    eaf = exposures[i, "eaf"],
    beta = exposures[i, "beta"],
    se = exposures[i, "se"],
    pval = exposures[i, "pval"],
    samplesize = exposures[i, "samplesize"],
    chr = exposures[i, "chr"],
    pos = exposures[i, "pos"],
    p_threshold = 5e-8,
    clump = TRUE,
    clump_kb = 10000,
    clump_r2 = 0.001
  )
}

# Extract outcome data ----
message("Extract outcome data")

outcomes <- gwas[gwas$type == "outcome", ]

for (i in 1:nrow(outcomes)) {
  out <- extract_outcome(
    instrument = exp$SNP,
    data = outcomes[i, "data"],
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
