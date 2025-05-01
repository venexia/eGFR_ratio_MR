# Specify gwas data frame ----
message("Specify gwas data frame")

gwas <- data.frame(
  type = character(),
  source = character(),
  data = character(),
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

# Add eGFR ----
message("Add eGFR")

gwas[nrow(gwas) + 1, ] <- c(
  type = "exposure",
  source = "https://ckdgen.imbi.uni-freiburg.de/files/Wuttke2019/20171017_MW_eGFR_overall_EA_nstud42.dbgap.txt.gz",
  data = "raw/20171017_MW_eGFR_overall_EA_nstud42.dbgap.txt.gz",
  phenotype = "eGFR",
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

# Add eGFR ratio ----
message("Add eGFR ratio")

gwas[nrow(gwas) + 1, ] <- c(
  type = "exposure",
  source = "https://ckdgen.imbi.uni-freiburg.de/files/Teumer2019/formatted_20180205-MA_overall-ALL-nstud_18-SumMac_400.tbl.rsid.gz",
  data = "raw/eGFR11_cys_cre_ratio_imputed.txt.gz",
  phenotype = "eGFR_ratio",
  snp = "SNP",
  effect_allele = "A1",
  other_allele = "ALLELE0",
  eaf = "A1FREQ",
  beta = "BETA",
  se = "SE",
  pval = "P_BOLT_LMM_INF",
  samplesize = "",
  chr = "CHR",
  pos = "BP"
)

# Add microalbuminuria ----
message("Add microalbuminuria")

gwas[nrow(gwas) + 1, ] <- c(
  type = "outcome",
  source = "https://ckdgen.imbi.uni-freiburg.de/files/Teumer2019/formatted_20180205-MA_overall-ALL-nstud_18-SumMac_400.tbl.rsid.gz",
  data = "raw/formatted_20180205-MA_overall-ALL-nstud_18-SumMac_400.tbl.rsid.gz",
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

gwas[nrow(gwas) + 1, ] <- c(
  type = "outcome",
  source = "https://ckdgen.imbi.uni-freiburg.de/files/Wuttke2019/CKD_overall_EA_JW_20180223_nstud23.dbgap.txt.gz",
  data = "raw/CKD_overall_EA_JW_20180223_nstud23.dbgap.txt.gz",
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
# gwas[nrow(gwas) + 1, ] <- c(
#   type = "outcome",
#   source = "",
#   data = "",
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
