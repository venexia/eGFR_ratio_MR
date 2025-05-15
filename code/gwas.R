# Specify gwas data frame ----
message("Specify gwas data frame")

gwas <- data.frame(
  type = character(),
  source = character(),
  data = character(),
  phenotype = character(),
  snp_col = character(),
  effect_allele_col = character(),
  other_allele_col = character(),
  eaf_col = character(),
  beta_col = character(),
  se_col = character(),
  pval_col = character(),
  samplesize_col = character(),
  samplesize = numeric(),
  chr_col = character(),
  pos_col = character(),
  build = character(),
  ukb = numeric()
)

# Add eGFR ratio ----
message("Add eGFR ratio")

gwas[nrow(gwas) + 1, ] <- c(
  type = "exposure",
  source = "MRC IEU pipeline",
  data = "raw/eGFR11_cys_cre_ratio_imputed.txt.gz",
  phenotype = "eGFR_ratio",
  snp = "SNP",
  effect_allele = "A1",
  other_allele = "ALLELE0",
  eaf = "A1FREQ",
  beta = "BETA",
  se = "SE",
  pval = "P_BOLT_LMM_INF",
  samplesize_col = "",
  samplesize = 502134, # Check with CB
  chr = "CHR",
  pos = "BP",
  build = "GRCh37",
  ukb = 1
)

# Add eGFR ----
message("Add eGFR")

gwas[nrow(gwas) + 1, ] <- c(
  type = "exposure",
  source = "https://ckdgen.imbi.uni-freiburg.de/files/Wuttke2019/20171017_MW_eGFR_overall_EA_nstud42.dbgap.txt.gz",
  data = "raw/20171017_MW_eGFR_overall_EA_nstud42.dbgap.txt.gz",
  phenotype = "eGFR",
  snp_col = "RSID",
  effect_allele_col = "Allele1",
  other_allele_col = "Allele2",
  eaf_col = "Freq1",
  beta_col = "Effect",
  se_col = "StdErr",
  pval_col = "P-value",
  samplesize_col = "n_total_sum",
  samplesize = 567460,
  chr_col = "Chr",
  pos_col = "Pos_b37",
  build = "GRCh37",
  ukb = 0 # UKB was not included but used to follow up SNPs
)

# Add microalbuminuria ----
message("Add microalbuminuria")

gwas[nrow(gwas) + 1, ] <- c(
  type = "outcome",
  source = "https://ckdgen.imbi.uni-freiburg.de/files/Teumer2019/formatted_20180205-MA_overall-ALL-nstud_18-SumMac_400.tbl.rsid.gz",
  data = "raw/formatted_20180205-MA_overall-ALL-nstud_18-SumMac_400.tbl.rsid.gz",
  phenotype = "Microalbuminuria",
  snp_col = "RSID",
  effect_allele_col = "Allele1",
  other_allele_col = "Allele2",
  eaf_col = "Freq1",
  beta_col = "Effect",
  se_col = "StdErr",
  pval_col = "P-value",
  samplesize_col = "n_total_sum",
  samplesize = 347269,
  chr_col = "Chr",
  pos_col = "Pos_b37",
  build = "GRCh37",
  ukb = (33296 / 347269)
)

# Add chronic kidney disease ----
message("Add chronic kidney disease")

gwas[nrow(gwas) + 1, ] <- c(
  type = "outcome",
  source = "https://ckdgen.imbi.uni-freiburg.de/files/Wuttke2019/CKD_overall_EA_JW_20180223_nstud23.dbgap.txt.gz",
  data = "raw/CKD_overall_EA_JW_20180223_nstud23.dbgap.txt.gz",
  phenotype = "Chronic kidney disease",
  snp_col = "RSID",
  effect_allele_col = "Allele1",
  other_allele_col = "Allele2",
  eaf_col = "Freq1",
  beta_col = "Effect",
  se_col = "StdErr",
  pval_col = "P-value",
  samplesize_col = "n_total_sum",
  samplesize = 480698,
  chr_col = "Chr",
  pos_col = "Pos_b37",
  build = "GRCh37",
  ukb = 0 # UKB was not included but used to follow up SNPs
)

# Add cardiovascular disease ----
message("Add cardiovascular disease")

gwas[nrow(gwas) + 1, ] <- c(
  type = "outcome",
  source = "https://personal.broadinstitute.org/ryank/Aragam_2022_CARDIoGRAM_CAD_GWAS.zip",
  data = "raw/CAD_GWAS_primary_discovery_meta.tsv",
  phenotype = "Cardiovascular disease",
  snp_col = "MarkerName",
  effect_allele_col = "Allele1",
  other_allele_col = "Allele2",
  eaf_col = "Freq1",
  beta = "Effect",
  se_col = "StdErr",
  pval_col = "P-value",
  samplesize_col = "N",
  samplesize = 1165690,
  chr_col = "CHR",
  pos_col = "BP",
  build = "GRCh37",
  ukb = (472335 / 1165690)
)
