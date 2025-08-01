# Load IEU Open GWAS catalog ----
message("Load IEU Open GWAS catalog")

ao <- TwoSampleMR::available_outcomes()

# Specify gwas data frame ----
message("Specify gwas data frame")

gwas <- data.frame(
  map_snps = character(),
  source = character(),
  data = character(),
  phenotype = character(),
  phenotype_short = character(),
  category = character(),
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
  adj = character(),
  ukb = numeric(),
  id = character()
)

# Add eGFR ratio ----
message("Add eGFR ratio")

gwas[nrow(gwas) + 1, ] <- c(
  map_snps = FALSE,
  source = "MRC IEU pipeline",
  data = "raw/eGFR11_cys_cre_ratio_imputed.txt.gz",
  phenotype = "eGFR_ratio",
  phenotype_short = "eGFRr",
  category = "continuous",
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
  adj = "",
  ukb = 1,
  id = ""
)

# Add eGFR ----
message("Add eGFR")

gwas[nrow(gwas) + 1, ] <- c(
  map_snps = FALSE,
  source = "https://ckdgen.imbi.uni-freiburg.de/files/Wuttke2019/20171017_MW_eGFR_overall_EA_nstud42.dbgap.txt.gz",
  data = "raw/20171017_MW_eGFR_overall_EA_nstud42.dbgap.txt.gz",
  phenotype = "eGFR",
  phenotype_short = "eGFR",
  category = "continuous",
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
  adj = "",
  ukb = 0, # UKB was not included but used to follow up SNPs
  id = ""
)

# Add microalbuminuria ----
message("Add microalbuminuria")

gwas[nrow(gwas) + 1, ] <- c(
  map_snps = FALSE,
  source = "https://ckdgen.imbi.uni-freiburg.de/files/Teumer2019/formatted_20180205-MA_overall-ALL-nstud_18-SumMac_400.tbl.rsid.gz",
  data = "raw/formatted_20180205-MA_overall-ALL-nstud_18-SumMac_400.tbl.rsid.gz",
  phenotype = "Microalbuminuria",
  phenotype_short = "MA",
  category = "binary",
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
  adj = "",
  ukb = (33296 / 347269),
  id = ""
)

# Add chronic kidney disease ----
message("Add chronic kidney disease")

gwas[nrow(gwas) + 1, ] <- c(
  map_snps = FALSE,
  source = "https://ckdgen.imbi.uni-freiburg.de/files/Wuttke2019/CKD_overall_EA_JW_20180223_nstud23.dbgap.txt.gz",
  data = "raw/CKD_overall_EA_JW_20180223_nstud23.dbgap.txt.gz",
  phenotype = "Chronic kidney disease",
  phenotype_short = "CKD",
  category = "binary",
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
  adj = "",
  ukb = 0, # UKB was not included but used to follow up SNPs
  id = ""
)

# Add cardiovascular disease ----
message("Add cardiovascular disease")

id <- "ebi-a-GCST005195"

gwas[nrow(gwas) + 1, ] <- c(
  map_snps = FALSE,
  source = "IEU Open GWAS",
  data = "",
  phenotype = "Coronary artery disease",
  phenotype_short = "CAD",
  category = "binary",
  snp_col = "",
  effect_allele_col = "",
  other_allele_col = "",
  eaf_col = "",
  beta = "",
  se_col = "",
  pval_col = "",
  samplesize_col = "",
  samplesize = ao[ao$id == id, ]$sample_size,
  chr_col = "",
  pos_col = "",
  build = "GRCh37",
  adj = "",
  ukb = (472335 / 1165690),
  id = id
)


# Add body mass index ----
message("Add body mass index")

id <- "ieu-b-40"

gwas[nrow(gwas) + 1, ] <- c(
  map_snps = FALSE,
  source = "IEU Open GWAS",
  data = "",
  phenotype = "Body mass index",
  phenotype_short = "BMI",
  category = "continuous",
  snp_col = "",
  effect_allele_col = "",
  other_allele_col = "",
  eaf_col = "",
  beta = "",
  se_col = "",
  pval_col = "",
  samplesize_col = "",
  samplesize = ao[ao$id == id, ]$sample_size,
  chr_col = "",
  pos_col = "",
  build = "GRCh37",
  adj = "",
  ukb = (456426 / ao[ao$id == id, ]$sample_size),
  id = id
)

# Add hip circumference ----
message("Add hip circumference")

id <- "ukb-b-15590"

gwas[nrow(gwas) + 1, ] <- c(
  map_snps = FALSE,
  source = "IEU Open GWAS",
  data = "",
  phenotype = "Hip circumference",
  phenotype_short = "HC",
  category = "continuous",
  snp_col = "",
  effect_allele_col = "",
  other_allele_col = "",
  eaf_col = "",
  beta = "",
  se_col = "",
  pval_col = "",
  samplesize_col = "",
  samplesize = ao[ao$id == id, ]$sample_size,
  chr_col = "",
  pos_col = "",
  build = "GRCh37",
  adj = "",
  ukb = 1,
  id = id
)

# Add waist circumference ----
message("Add waist circumference")

id <- "ukb-b-9405"

gwas[nrow(gwas) + 1, ] <- c(
  map_snps = FALSE,
  source = "IEU Open GWAS",
  data = "",
  phenotype = "Waist circumference",
  phenotype_short = "WC",
  category = "continuous",
  snp_col = "",
  effect_allele_col = "",
  other_allele_col = "",
  eaf_col = "",
  beta = "",
  se_col = "",
  pval_col = "",
  samplesize_col = "",
  samplesize = ao[ao$id == id, ]$sample_size,
  chr_col = "",
  pos_col = "",
  build = "GRCh37",
  adj = "",
  ukb = 1,
  id = id
)

# Add HbA1c ----
message("Add HbA1c")

id <- "ieu-b-104"

gwas[nrow(gwas) + 1, ] <- c(
  map_snps = FALSE,
  source = "IEU Open GWAS",
  data = "",
  phenotype = "HbA1c",
  phenotype_short = "HbA1c",
  category = "continuous",
  snp_col = "",
  effect_allele_col = "",
  other_allele_col = "",
  eaf_col = "",
  beta = "",
  se_col = "",
  pval_col = "",
  samplesize_col = "",
  samplesize = ao[ao$id == id, ]$sample_size,
  chr_col = "",
  pos_col = "",
  build = "GRCh37",
  adj = "",
  ukb = 0,
  id = id
)

# Add diastolic blood pressure ----
message("Add diastolic blood pressure")

id <- "ebi-a-GCST90000063"

gwas[nrow(gwas) + 1, ] <- c(
  map_snps = FALSE,
  source = "IEU Open GWAS",
  data = "",
  phenotype = "Diastolic blood pressure",
  phenotype_short = "DBP",
  category = "continuous",
  snp_col = "",
  effect_allele_col = "",
  other_allele_col = "",
  eaf_col = "",
  beta = "",
  se_col = "",
  pval_col = "",
  samplesize_col = "",
  samplesize = ao[ao$id == id, ]$sample_size,
  chr_col = "",
  pos_col = "",
  build = "GRCh37",
  adj = "",
  ukb = (364510 / ao[ao$id == id, ]$sample_size),
  id = id
)

# Add systolic blood pressure ----
message("Add systolic blood pressure")

id <- "ebi-a-GCST90000062"

gwas[nrow(gwas) + 1, ] <- c(
  map_snps = FALSE,
  source = "IEU Open GWAS",
  data = "",
  phenotype = "Systolic blood pressure",
  phenotype_short = "SBP",
  category = "continuous",
  snp_col = "",
  effect_allele_col = "",
  other_allele_col = "",
  eaf_col = "",
  beta = "",
  se_col = "",
  pval_col = "",
  samplesize_col = "",
  samplesize = ao[ao$id == id, ]$sample_size,
  chr_col = "",
  pos_col = "",
  build = "GRCh37",
  adj = "",
  ukb = (364510 / ao[ao$id == id, ]$sample_size),
  id = id
)

# Add smoking ----
message("Add smoking")

gwas[nrow(gwas) + 1, ] <- c(
  map_snps = FALSE,
  source = "https://doi.org/10.5523/bris.10i96zb8gm0j81yz0q6ztei23d",
  data = "raw/2019.10.02 Lifetime Smoking GWAS Data Sheet 1.txt",
  phenotype = "Lifetime smoking index",
  phenotype_short = "LSI",
  category = "continuous",
  snp_col = "SNP",
  effect_allele_col = "EFFECT_ALLELE",
  other_allele_col = "OTHER_ALLELE",
  eaf_col = "EAF",
  beta_col = "BETA",
  se_col = "SE",
  pval_col = "P",
  samplesize_col = "",
  samplesize = 462690,
  chr_col = "CHR",
  pos_col = "BP",
  build = "GRCh37",
  adj = "",
  ukb = 1, # UKB was not included but used to follow up SNPs
  id = ""
)

# Add alcohol ----
message("Add alcohol")

id <- "ieu-a-1283"

gwas[nrow(gwas) + 1, ] <- c(
  map_snps = FALSE,
  source = "IEU Open GWAS",
  data = "",
  phenotype = "Alcohol consumption",
  phenotype_short = "Alc",
  category = "continuous",
  snp_col = "",
  effect_allele_col = "",
  other_allele_col = "",
  eaf_col = "",
  beta = "",
  se_col = "",
  pval_col = "",
  samplesize_col = "",
  samplesize = ao[ao$id == id, ]$sample_size,
  chr_col = "",
  pos_col = "",
  build = "GRCh37",
  adj = "",
  ukb = 1,
  id = id
)

# Add hypertension ----
message("Add hypertension")

id <- "ukb-d-I9_HYPTENS"

gwas[nrow(gwas) + 1, ] <- c(
  map_snps = FALSE,
  source = "IEU Open GWAS",
  data = "",
  phenotype = "Hypertension",
  phenotype_short = "Hypertension",
  category = "binary",
  snp_col = "",
  effect_allele_col = "",
  other_allele_col = "",
  eaf_col = "",
  beta = "",
  se_col = "",
  pval_col = "",
  samplesize_col = "",
  samplesize = ao[ao$id == id, ]$sample_size,
  chr_col = "",
  pos_col = "",
  build = "GRCh37",
  adj = "",
  ukb = 1,
  id = id
)

# # Add hypotension ----
# message("Add hypotension")
#
# id <- "finn-b-I9_HYPOTE"
#
# gwas[nrow(gwas) + 1, ] <- c(
#   map_snps = FALSE,
#   source = "IEU Open GWAS",
#   data = "",
#   phenotype = "Hypotension",
#   phenotype_short = "Hypotension",
#   category = "binary",
#   snp_col = "",
#   effect_allele_col = "",
#   other_allele_col = "",
#   eaf_col = "",
#   beta = "",
#   se_col = "",
#   pval_col = "",
#   samplesize_col = "",
#   samplesize = (1857 + 216463),
#   chr_col = "",
#   pos_col = "",
#   build = "GRCh37",
#   adj = "",
#   ukb = 0,
#   id = id
# )

# Add atrial fibrillation ----
message("Add atrial fibrillation")

id <- "ebi-a-GCST006414"

gwas[nrow(gwas) + 1, ] <- c(
  map_snps = FALSE,
  source = "IEU Open GWAS",
  data = "",
  phenotype = "Atrial fibrillation",
  category = "binary",
  phenotype_short = "AF",
  snp_col = "",
  effect_allele_col = "",
  other_allele_col = "",
  eaf_col = "",
  beta = "",
  se_col = "",
  pval_col = "",
  samplesize_col = "",
  samplesize = ao[ao$id == id, ]$sample_size,
  chr_col = "",
  pos_col = "",
  build = "GRCh37",
  adj = "",
  ukb = (395739 / 1030836),
  id = id
)

# Add chronic ischaemic heart disease ----
message("Add chronic ischaemic heart disease")

id <- "ukb-d-I9_IHD"

gwas[nrow(gwas) + 1, ] <- c(
  map_snps = FALSE,
  source = "IEU Open GWAS",
  data = "",
  phenotype = "Ischaemic heart disease",
  phenotype_short = "IHD",
  category = "binary",
  snp_col = "",
  effect_allele_col = "",
  other_allele_col = "",
  eaf_col = "",
  beta = "",
  se_col = "",
  pval_col = "",
  samplesize_col = "",
  samplesize = ao[ao$id == id, ]$sample_size,
  chr_col = "",
  pos_col = "",
  build = "GRCh37",
  adj = "",
  ukb = 1,
  id = id
)

# # Add type 2 diabetes ----
# message("Add type 2 diabetes")
#
# gwas[nrow(gwas) + 1, ] <- c(
#   map_snps = FALSE,
#   source = "https://diagram-consortium.org/downloads.html; 'T2D GWAS meta-analysis - Summary of T2D associations, unadjusted for BMI and without UK Biobank subjects'",
#   data = "raw/Mahajan.NatGenet2018b.T2D-noUKBB.European.zip",
#   phenotype = "Type 2 diabetes",
#   phenotype_short = "T2D",
#   category = "binary",
#   snp_col = "SNP",
#   effect_allele_col = "EA",
#   other_allele_col = "NEA",
#   eaf_col = "EAF",
#   beta_col = "Beta",
#   se_col = "SE",
#   pval_col = "Pvalue",
#   samplesize_col = "",
#   samplesize = 898130,
#   chr_col = "Chr",
#   pos_col = "Pos",
#   build = "GRCh37",
#   adj = "",
#   ukb = 0,
#   id = ""
# )

# Add renal failure ----
message("Add renal failure")

id <- "ebi-a-GCST90018822"

gwas[nrow(gwas) + 1, ] <- c(
  map_snps = FALSE,
  source = "IEU Open GWAS",
  data = "",
  phenotype = "Chronic renal failure",
  phenotype_short = "RF",
  category = "binary",
  snp_col = "",
  effect_allele_col = "",
  other_allele_col = "",
  eaf_col = "",
  beta = "",
  se_col = "",
  pval_col = "",
  samplesize_col = "",
  samplesize = ao[ao$id == id, ]$sample_size,
  chr_col = "",
  pos_col = "",
  build = "GRCh37",
  adj = "",
  ukb = ((6542 + 342005) / ao[ao$id == id, ]$sample_size),
  id = id
)

# Add cerebral infarction ----
message("Add cerebral infarction")

id <- "ukb-d-I63"

gwas[nrow(gwas) + 1, ] <- c(
  map_snps = FALSE,
  source = "IEU Open GWAS",
  data = "",
  phenotype = "Cerebral infarction",
  phenotype_short = "CI",
  category = "binary",
  snp_col = "",
  effect_allele_col = "",
  other_allele_col = "",
  eaf_col = "",
  beta = "",
  se_col = "",
  pval_col = "",
  samplesize_col = "",
  samplesize = ao[ao$id == id, ]$sample_size,
  chr_col = "",
  pos_col = "",
  build = "GRCh37",
  adj = "",
  ukb = 1,
  id = id
)

# Add type 1 diabetes ----
message("Add type 1 diabetes")

id <- "ebi-a-GCST90014023"

gwas[nrow(gwas) + 1, ] <- c(
  map_snps = FALSE,
  source = "IEU Open GWAS",
  data = "",
  phenotype = "Type 1 diabetes",
  phenotype_short = "T1D",
  category = "binary",
  snp_col = "",
  effect_allele_col = "",
  other_allele_col = "",
  eaf_col = "",
  beta = "",
  se_col = "",
  pval_col = "",
  samplesize_col = "",
  samplesize = ao[ao$id == id, ]$sample_size,
  chr_col = "",
  pos_col = "",
  build = "GRCh37",
  adj = "",
  ukb = (363495 / ao[ao$id == id, ]$sample_size),
  id = id
)

# Add angina pectoris ----
message("Add angina pectoris")

id <- "ebi-a-GCST90018793"

gwas[nrow(gwas) + 1, ] <- c(
  map_snps = FALSE,
  source = "IEU Open GWAS",
  data = "",
  phenotype = "Angina pectoris",
  phenotype_short = "AP",
  category = "binary",
  snp_col = "",
  effect_allele_col = "",
  other_allele_col = "",
  eaf_col = "",
  beta = "",
  se_col = "",
  pval_col = "",
  samplesize_col = "",
  samplesize = ao[ao$id == id, ]$sample_size,
  chr_col = "",
  pos_col = "",
  build = "GRCh37",
  adj = "",
  ukb = ((17894 + 325132) / ao[ao$id == id, ]$sample_size),
  id = id
)

# Add acute myocardial infarction ----
message("Add acute myocardial infarction")

id <- "ebi-a-GCST90018877"

gwas[nrow(gwas) + 1, ] <- c(
  map_snps = FALSE,
  source = "IEU Open GWAS",
  data = "",
  phenotype = "Myocardial infarction",
  phenotype_short = "MI",
  category = "binary",
  snp_col = "",
  effect_allele_col = "",
  other_allele_col = "",
  eaf_col = "",
  beta = "",
  se_col = "",
  pval_col = "",
  samplesize_col = "",
  samplesize = ao[ao$id == id, ]$sample_size,
  chr_col = "",
  pos_col = "",
  build = "GRCh37",
  adj = "",
  ukb = ((8234 + 115774) / ao[ao$id == id, ]$sample_size),
  id = id
)

# Add heart failure ----
message("Add heart failure")

id <- "ebi-a-GCST009541"

gwas[nrow(gwas) + 1, ] <- c(
  map_snps = FALSE,
  source = "IEU Open GWAS",
  data = "",
  phenotype = "Heart failure",
  phenotype_short = "HF",
  category = "binary",
  snp_col = "",
  effect_allele_col = "",
  other_allele_col = "",
  eaf_col = "",
  beta = "",
  se_col = "",
  pval_col = "",
  samplesize_col = "",
  samplesize = ao[ao$id == id, ]$sample_size,
  chr_col = "",
  pos_col = "",
  build = "GRCh37",
  adj = "",
  ukb = (488377 / ao[ao$id == id, ]$sample_size),
  id = id
)

# Add atherosclerosis ----
message("Add atherosclerosis")

id <- "ukb-d-I9_CORATHER"

gwas[nrow(gwas) + 1, ] <- c(
  map_snps = FALSE,
  source = "IEU Open GWAS",
  data = "",
  phenotype = "Coronary atherosclerosis",
  phenotype_short = "CA",
  category = "binary",
  snp_col = "",
  effect_allele_col = "",
  other_allele_col = "",
  eaf_col = "",
  beta = "",
  se_col = "",
  pval_col = "",
  samplesize_col = "",
  samplesize = ao[ao$id == id, ]$sample_size,
  chr_col = "",
  pos_col = "",
  build = "GRCh37",
  adj = "",
  ukb = 1,
  id = id
)

# Add aortic aneurysm ----
message("Add aortic aneurysm")

id <- "ebi-a-GCST90018783"

gwas[nrow(gwas) + 1, ] <- c(
  map_snps = FALSE,
  source = "IEU Open GWAS",
  data = "",
  phenotype = "Aortic aneurysm",
  phenotype_short = "AA",
  category = "binary",
  snp_col = "",
  effect_allele_col = "",
  other_allele_col = "",
  eaf_col = "",
  beta = "",
  se_col = "",
  pval_col = "",
  samplesize_col = "",
  samplesize = ao[ao$id == id, ]$sample_size,
  chr_col = "",
  pos_col = "",
  build = "GRCh37",
  adj = "",
  ukb = ((1520 + 128425) / ao[ao$id == id, ]$sample_size),
  id = id
)

# # Add arterial embolism and thrombosis ----
# message("Add arterial embolism and thrombosis")
#
# id <- "finn-b-I9_ARTEMBTHR"
#
# gwas[nrow(gwas) + 1, ] <- c(
#   map_snps = FALSE,
#   source = "IEU Open GWAS",
#   data = "",
#   phenotype = "Arterial embolism and thrombosis",
#   phenotype_short = "AET",
#   category = "binary",
#   snp_col = "",
#   effect_allele_col = "",
#   other_allele_col = "",
#   eaf_col = "",
#   beta = "",
#   se_col = "",
#   pval_col = "",
#   samplesize_col = "",
#   samplesize = (789 + 206541),
#   chr_col = "",
#   pos_col = "",
#   build = "GRCh37",
#   adj = "",
#   ukb = 0,
#   id = id
# )

# Add pulmonary embolism ----
message("Add pulmonary embolism")

id <- "finn-b-I9_PULMEMB"

gwas[nrow(gwas) + 1, ] <- c(
  map_snps = FALSE,
  source = "IEU Open GWAS",
  data = "",
  phenotype = "Pulmonary embolism",
  phenotype_short = "PE",
  category = "binary",
  snp_col = "",
  effect_allele_col = "",
  other_allele_col = "",
  eaf_col = "",
  beta = "",
  se_col = "",
  pval_col = "",
  samplesize_col = "",
  samplesize = (4185 + 214228),
  chr_col = "",
  pos_col = "",
  build = "GRCh37",
  adj = "",
  ukb = 0,
  id = id
)

# Format ----

gwas$samplesize <- as.numeric(gwas$samplesize)
gwas$ukb <- as.numeric(gwas$ukb)

# Save ----

data.table::fwrite(gwas, "data/gwas.csv")
