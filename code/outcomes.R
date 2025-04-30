# Specify outcomes data frame ----
message("Specify outcomes data frame")

outcomes <- data.frame(
  source = character(),
  path = character(),
  phenotype = character(),
  snp = character(),
  effect_allele = character(),
  other_allele = character(),
  eaf = character(),
  beta = character(),
  se = character(),
  pval = character(),
  samplesize = character(),
  chr = character(),
  pos = character()
)

# Add microalbuminuria ----
message("Add microalbuminuria")

outcomes[nrow(outcomes) + 1, ] <- c(
  source = "https://ckdgen.imbi.uni-freiburg.de/files/Teumer2019/formatted_20180205-MA_overall-ALL-nstud_18-SumMac_400.tbl.rsid.gz",
  path = "raw/formatted_20180205-MA_overall-ALL-nstud_18-SumMac_400.tbl.rsid.gz",
  phenotype = "Microalbuminuria",
  snp = "RSID",
  effect_allele = "Allele1",
  other_allele = "Allele2",
  eaf = "Freq1",
  beta = "Effect",
  se = "StdErr",
  pval = "P-value",
  samplesize = "n_total_sum",
  chr = "Chr",
  pos = "Pos_b37"
)

# Add chronic kidney disease ----
message("Add chronic kidney disease")

outcomes[nrow(outcomes) + 1, ] <- c(
  source = "https://ckdgen.imbi.uni-freiburg.de/files/Wuttke2019/CKD_overall_EA_JW_20180223_nstud23.dbgap.txt.gz",
  path = "raw/CKD_overall_EA_JW_20180223_nstud23.dbgap.txt.gz",
  phenotype = "Chronic kidney disease",
  snp = "RSID",
  effect_allele = "Allele1",
  other_allele = "Allele2",
  eaf = "Freq1",
  beta = "Effect",
  se = "StdErr",
  pval = "P-value",
  samplesize = "n_total_sum",
  chr = "Chr",
  pos = "Pos_b37"
)

# Add cardiovascular disease ----
# message("Add cardiovascular disease")
#
# outcomes[nrow(outcomes) + 1, ] <- c(
#   source = "",
#   path = "",
#   phenotype = "Cardiovascular disease",
#   snp = "RSID",
#   effect_allele = "Allele1",
#   other_allele = "Allele2",
#   eaf = "Freq1",
#   beta = "Effect",
#   se = "StdErr",
#   pval = "P-value",
#   samplesize = "n_total_sum",
#   chr = "Chr",
#   pos = "Pos_b37"
# )
