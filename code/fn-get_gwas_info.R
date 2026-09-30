get_gwas_info <- function() {
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
    ncase_col = character(),
    ncase = numeric(),
    ncontrol_col = character(),
    ncontrol = numeric(),
    chr_col = character(),
    pos_col = character(),
    build = character(),
    adj = character(),
    ukb = numeric(),
    id = character(),
    unit = character(),
    SD = numeric()
  )

  # Add eGFR ratio ----
  message("Add eGFR ratio")

  gwas[nrow(gwas) + 1, ] <- c(
    map_snps = "",
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
    samplesize = 502134, # Check with Chris or Sam!
    ncase_col = "",
    ncase = NA,
    ncontrol_col = "",
    ncontrol = NA,
    chr = "CHR",
    pos = "BP",
    build = "GRCh37",
    adj = "",
    ukb = 1,
    id = "",
    unit = "Unitless [Ratio]",
    SD = NA
  )

  # Add eGFR ----
  message("Add eGFR")

  gwas[nrow(gwas) + 1, ] <- c(
    map_snps = "",
    source = "https://ckdgen.imbi.uni-freiburg.de/files/Wuttke2019/20171017_MW_eGFR_overall_EA_nstud42.dbgap.txt.gz",
    data = "raw/20171017_MW_eGFR_overall_EA_nstud42.dbgap.txt.gz",
    phenotype = "eGFR (creatinine)",
    phenotype_short = "eGFRcrea",
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
    ncase_col = "",
    ncase = NA,
    ncontrol_col = "",
    ncontrol = NA,
    chr_col = "Chr",
    pos_col = "Pos_b37",
    build = "GRCh37",
    adj = "",
    ukb = 0, # UKB was not included but used to follow up SNPs
    id = "",
    unit = "mL/min/1.73 m^2",
    SD = 9.6 # Back-calculated from median 89 ml min−1 per 1.73 m2 (interquartile range, IQR: 81, 94)
  )

  gwas[nrow(gwas) + 1, ] <- c(
    map_snps = "",
    source = "https://ckdgen.imbi.uni-freiburg.de/files/Stanzick2021/metal_eGFRcys_meta1.TBL.map.annot.gc.gz",
    data = "raw/metal_eGFRcys_meta1.TBL.map.annot.gc.gz",
    phenotype = "eGFR (cystatin C)",
    phenotype_short = "eGFRcys",
    category = "continuous",
    snp_col = "RSID",
    effect_allele_col = "Allele1",
    other_allele_col = "Allele2",
    eaf_col = "Freq1",
    beta_col = "Effect",
    se_col = "StdErr",
    pval_col = "P.value",
    samplesize_col = "n",
    samplesize = 460826,
    ncase_col = "",
    ncase = NA,
    ncontrol_col = "",
    ncontrol = NA,
    chr_col = "chr",
    pos_col = "pos",
    build = "GRCh37",
    adj = "",
    ukb = 1,
    id = "",
    unit = "mL/min/1.73 m^2",
    SD = 18.7 # Reported in paper for HUNT
  )

  # Add microalbuminuria ----
  message("Add microalbuminuria")

  gwas[nrow(gwas) + 1, ] <- c(
    map_snps = "",
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
    ncase_col = "",
    ncase = 51861,
    ncontrol_col = "",
    ncontrol = 297093,
    chr_col = "Chr",
    pos_col = "Pos_b37",
    build = "GRCh37",
    adj = "",
    ukb = (33296 / 347269),
    id = "",
    unit = "Binary",
    SD = NA
  )

  # Add chronic kidney disease ----
  message("Add chronic kidney disease")

  gwas[nrow(gwas) + 1, ] <- c(
    map_snps = "",
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
    ncase_col = "",
    ncase = 41395,
    ncontrol_col = "",
    ncontrol = 439303,
    chr_col = "Chr",
    pos_col = "Pos_b37",
    build = "GRCh37",
    adj = "",
    ukb = 0, # UKB was not included but used to follow up SNPs
    id = "",
    unit = "Binary",
    SD = NA
  )

  # Add cardiovascular disease ----
  message("Add cardiovascular disease")

  id <- "ebi-a-GCST005195"

  gwas[nrow(gwas) + 1, ] <- c(
    map_snps = "",
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
    ncase_col = "",
    ncase = ao[ao$id == id, ]$sample_size - ao[ao$id == id, ]$ncontrol,
    ncontrol_col = "",
    ncontrol = ao[ao$id == id, ]$ncontrol,
    chr_col = "",
    pos_col = "",
    build = "GRCh37",
    adj = "",
    ukb = (472335 / 1165690),
    id = id,
    unit = "Binary",
    SD = NA
  )

  # Add body mass index ----
  message("Add body mass index")

  id <- "ieu-b-40"

  gwas[nrow(gwas) + 1, ] <- c(
    map_snps = "",
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
    ncase_col = "",
    ncase = NA,
    ncontrol_col = "",
    ncontrol = NA,
    chr_col = "",
    pos_col = "",
    build = "GRCh37",
    adj = "",
    ukb = (456426 / ao[ao$id == id, ]$sample_size),
    id = id,
    unit = "kg / m^2",
    SD = 4.76116 # From UK Biobank showcase (Data-Field 21001)
  )

  # Add hip circumference ----
  message("Add hip circumference")

  id <- "ukb-b-15590"

  gwas[nrow(gwas) + 1, ] <- c(
    map_snps = "",
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
    ncase_col = "",
    ncase = NA,
    ncontrol_col = "",
    ncontrol = NA,
    chr_col = "",
    pos_col = "",
    build = "GRCh37",
    adj = "",
    ukb = 1,
    id = id,
    unit = "SD",
    SD = NA
  )

  # Add waist circumference ----
  message("Add waist circumference")

  id <- "ukb-b-9405"

  gwas[nrow(gwas) + 1, ] <- c(
    map_snps = "",
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
    ncase_col = "",
    ncase = NA,
    ncontrol_col = "",
    ncontrol = NA,
    chr_col = "",
    pos_col = "",
    build = "GRCh37",
    adj = "",
    ukb = 1,
    id = id,
    unit = "SD",
    SD = NA
  )

  # Add whole body fat free mass ----
  message("Add whole body fat free mass")

  id <- "ukb-b-13354"

  gwas[nrow(gwas) + 1, ] <- c(
    map_snps = "",
    source = "IEU Open GWAS",
    data = "",
    phenotype = "Whole body fat-free mass",
    phenotype_short = "WBFFM",
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
    ncase_col = "",
    ncase = NA,
    ncontrol_col = "",
    ncontrol = NA,
    chr_col = "",
    pos_col = "",
    build = "GRCh37",
    adj = "",
    ukb = 1,
    id = id,
    unit = "SD",
    SD = NA
  )

  # Add leg fat free mass (right) ----
  message("Add leg fat free mass (right)")

  id <- "ukb-b-12828"

  gwas[nrow(gwas) + 1, ] <- c(
    map_snps = "",
    source = "IEU Open GWAS",
    data = "",
    phenotype = "Leg fat-free mass (right)",
    phenotype_short = "LFFMR",
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
    ncase_col = "",
    ncase = NA,
    ncontrol_col = "",
    ncontrol = NA,
    chr_col = "",
    pos_col = "",
    build = "GRCh37",
    adj = "",
    ukb = 1,
    id = id,
    unit = "SD",
    SD = NA
  )

  # Add HbA1c ----
  message("Add HbA1c")

  id <- "ieu-b-104"

  gwas[nrow(gwas) + 1, ] <- c(
    map_snps = "",
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
    ncase_col = "",
    ncase = NA,
    ncontrol_col = "",
    ncontrol = NA,
    chr_col = "",
    pos_col = "",
    build = "GRCh37",
    adj = "",
    ukb = 0,
    id = id,
    unit = "Percentage",
    SD = ao[ao$id == id, ]$sd
  )

  # Add diastolic blood pressure ----
  message("Add diastolic blood pressure")

  id <- "ebi-a-GCST90000063"

  gwas[nrow(gwas) + 1, ] <- c(
    map_snps = "",
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
    ncase_col = "",
    ncase = NA,
    ncontrol_col = "",
    ncontrol = NA,
    chr_col = "",
    pos_col = "",
    build = "GRCh37",
    adj = "age,age_sq,sex,bmi,medication,PCs",
    ukb = (364510 / ao[ao$id == id, ]$sample_size),
    id = id,
    unit = "Unitless [IRNT]",
    SD = NA
  )

  # Add systolic blood pressure ----
  message("Add systolic blood pressure")

  id <- "ebi-a-GCST90000062"

  gwas[nrow(gwas) + 1, ] <- c(
    map_snps = "",
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
    ncase_col = "",
    ncase = NA,
    ncontrol_col = "",
    ncontrol = NA,
    chr_col = "",
    pos_col = "",
    build = "GRCh37",
    adj = "age,age_sq,sex,bmi,medication,PCs",
    ukb = (364510 / ao[ao$id == id, ]$sample_size),
    id = id,
    unit = "Unitless [IRNT]",
    SD = NA
  )

  # Add smoking ----
  message("Add smoking")

  gwas[nrow(gwas) + 1, ] <- c(
    map_snps = "",
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
    ncase_col = "",
    ncase = NA,
    ncontrol_col = "",
    ncontrol = NA,
    chr_col = "CHR",
    pos_col = "BP",
    build = "GRCh37",
    adj = "",
    ukb = 1, # UKB was not included but used to follow up SNPs
    id = "",
    unit = "Per score unit",
    SD = 0.6940093
  )

  # Add alcohol ----
  message("Add alcohol")

  id <- "ieu-a-1283"

  gwas[nrow(gwas) + 1, ] <- c(
    map_snps = "",
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
    ncase_col = "",
    ncase = NA,
    ncontrol_col = "",
    ncontrol = NA,
    chr_col = "",
    pos_col = "",
    build = "GRCh37",
    adj = "",
    ukb = 1,
    id = id,
    unit = "Unitless [Residual log(alcohol units + 1)]",
    SD = NA
  )

  # Add hypertension ----
  message("Add hypertension")

  id <- "ieu-b-5144"

  gwas[nrow(gwas) + 1, ] <- c(
    map_snps = "",
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
    ncase_col = "",
    ncase = ao[ao$id == id, ]$sample_size - ao[ao$id == id, ]$ncontrol,
    ncontrol_col = "",
    ncontrol = ao[ao$id == id, ]$ncontrol,
    chr_col = "",
    pos_col = "",
    build = "GRCh37",
    adj = "",
    ukb = 1,
    id = id,
    unit = "Binary",
    SD = NA
  )

  # # Add hypotension ----
  # message("Add hypotension")
  #
  # id <- "finn-b-I9_HYPOTE"
  #
  # gwas[nrow(gwas) + 1, ] <- c(
  #   map_snps = "",
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
  #   ncase_col = "",
  #   ncase = ao[ao$id == id, ]$sample_size - ao[ao$id == id, ]$ncontrol,
  #   ncontrol_col = "",
  #   ncontrol = ao[ao$id == id, ]$ncontrol,
  #   chr_col = "",
  #   pos_col = "",
  #   build = "GRCh37",
  #   adj = "",
  #   ukb = 0,
  #   id = id,
  #   unit = "Binary",
  #   SD = NA
  # )

  # Add atrial fibrillation ----
  message("Add atrial fibrillation")

  id <- "ebi-a-GCST006414"

  gwas[nrow(gwas) + 1, ] <- c(
    map_snps = "",
    source = "IEU Open GWAS",
    data = "",
    phenotype = "Atrial fibrillation",
    phenotype_short = "AF",
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
    ncase_col = "",
    ncase = ao[ao$id == id, ]$sample_size - ao[ao$id == id, ]$ncontrol,
    ncontrol_col = "",
    ncontrol = ao[ao$id == id, ]$ncontrol,
    chr_col = "",
    pos_col = "",
    build = "GRCh37",
    adj = "",
    ukb = (395739 / 1030836),
    id = id,
    unit = "Binary",
    SD = NA
  )

  # Add chronic ischaemic heart disease ----
  message("Add chronic ischaemic heart disease")

  id <- "ukb-d-I9_IHD"

  gwas[nrow(gwas) + 1, ] <- c(
    map_snps = "",
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
    ncase_col = "",
    ncase = ao[ao$id == id, ]$sample_size - ao[ao$id == id, ]$ncontrol,
    ncontrol_col = "",
    ncontrol = ao[ao$id == id, ]$ncontrol,
    chr_col = "",
    pos_col = "",
    build = "GRCh37",
    adj = "",
    ukb = 1,
    id = id,
    unit = "Binary",
    SD = NA
  )

  # Add type 2 diabetes ----
  message("Add type 2 diabetes")

  gwas[nrow(gwas) + 1, ] <- c(
    map_snps = "chr;pos;effect_allele;other_allele",
    source = "https://www.diagram-consortium.org/downloads.html; 'T2DGGI GWAS EUR ancestry meta-analysis summary statistics'",
    data = "raw/Suzuki.Nature2024.T2DGGI.EUR.sumstats.zip",
    phenotype = "Type 2 diabetes",
    phenotype_short = "T2D",
    category = "binary",
    snp_col = "", # Sort
    effect_allele_col = "EffectAllele",
    other_allele_col = "NonEffectAllele",
    eaf_col = "EAF",
    beta_col = "Beta",
    se_col = "SE",
    pval_col = "Pval",
    samplesize_col = "Neff",
    samplesize = (242283 + 1569734),
    ncase_col = "",
    ncase = 428452,
    ncontrol_col = "",
    ncontrol = 2107149,
    chr_col = "Chromsome",
    pos_col = "Position",
    build = "GRCh37",
    adj = "",
    ukb = 1,
    id = "",
    unit = "Binary",
    SD = NA
  )

  # Add renal failure ----
  message("Add renal failure")

  id <- "ebi-a-GCST90018822"

  gwas[nrow(gwas) + 1, ] <- c(
    map_snps = "",
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
    ncase_col = "",
    ncase = ao[ao$id == id, ]$sample_size - ao[ao$id == id, ]$ncontrol,
    ncontrol_col = "",
    ncontrol = ao[ao$id == id, ]$ncontrol,
    chr_col = "",
    pos_col = "",
    build = "GRCh37",
    adj = "",
    ukb = ((6542 + 342005) / ao[ao$id == id, ]$sample_size),
    id = id,
    unit = "Binary",
    SD = NA
  )

  # Add cerebral infarction ----
  message("Add cerebral infarction")

  id <- "ukb-d-I63"

  gwas[nrow(gwas) + 1, ] <- c(
    map_snps = "",
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
    ncase_col = "",
    ncase = ao[ao$id == id, ]$sample_size - ao[ao$id == id, ]$ncontrol,
    ncontrol_col = "",
    ncontrol = ao[ao$id == id, ]$ncontrol,
    chr_col = "",
    pos_col = "",
    build = "GRCh37",
    adj = "",
    ukb = 1,
    id = id,
    unit = "Binary",
    SD = NA
  )

  # Add type 1 diabetes ----
  message("Add type 1 diabetes")

  id <- "ebi-a-GCST90014023"

  gwas[nrow(gwas) + 1, ] <- c(
    map_snps = "",
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
    ncase_col = "",
    ncase = ao[ao$id == id, ]$sample_size - ao[ao$id == id, ]$ncontrol,
    ncontrol_col = "",
    ncontrol = ao[ao$id == id, ]$ncontrol,
    chr_col = "",
    pos_col = "",
    build = "GRCh37",
    adj = "",
    ukb = (363495 / ao[ao$id == id, ]$sample_size),
    id = id,
    unit = "Binary",
    SD = NA
  )

  # Add angina pectoris ----
  message("Add angina pectoris")

  id <- "ebi-a-GCST90018793"

  gwas[nrow(gwas) + 1, ] <- c(
    map_snps = "",
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
    ncase_col = "",
    ncase = ao[ao$id == id, ]$sample_size - ao[ao$id == id, ]$ncontrol,
    ncontrol_col = "",
    ncontrol = ao[ao$id == id, ]$ncontrol,
    chr_col = "",
    pos_col = "",
    build = "GRCh37",
    adj = "",
    ukb = ((17894 + 325132) / ao[ao$id == id, ]$sample_size),
    id = id,
    unit = "Binary",
    SD = NA
  )

  # Add acute myocardial infarction ----
  message("Add acute myocardial infarction")

  id <- "ebi-a-GCST90018877"

  gwas[nrow(gwas) + 1, ] <- c(
    map_snps = "",
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
    ncase_col = "",
    ncase = ao[ao$id == id, ]$sample_size - ao[ao$id == id, ]$ncontrol,
    ncontrol_col = "",
    ncontrol = ao[ao$id == id, ]$ncontrol,
    chr_col = "",
    pos_col = "",
    build = "GRCh37",
    adj = "",
    ukb = ((8234 + 115774) / ao[ao$id == id, ]$sample_size),
    id = id,
    unit = "Binary",
    SD = NA
  )

  # Add heart failure ----
  message("Add heart failure")

  id <- "ebi-a-GCST009541"

  gwas[nrow(gwas) + 1, ] <- c(
    map_snps = "",
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
    ncase_col = "",
    ncase = ao[ao$id == id, ]$sample_size - ao[ao$id == id, ]$ncontrol,
    ncontrol_col = "",
    ncontrol = ao[ao$id == id, ]$ncontrol,
    chr_col = "",
    pos_col = "",
    build = "GRCh37",
    adj = "",
    ukb = (488377 / ao[ao$id == id, ]$sample_size),
    id = id,
    unit = "Binary",
    SD = NA
  )

  # Add atherosclerosis ----
  message("Add atherosclerosis")

  id <- "ukb-d-I9_CORATHER"

  gwas[nrow(gwas) + 1, ] <- c(
    map_snps = "",
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
    ncase_col = "",
    ncase = ao[ao$id == id, ]$sample_size - ao[ao$id == id, ]$ncontrol,
    ncontrol_col = "",
    ncontrol = ao[ao$id == id, ]$ncontrol,
    chr_col = "",
    pos_col = "",
    build = "GRCh37",
    adj = "",
    ukb = 1,
    id = id,
    unit = "Binary",
    SD = NA
  )

  # Add aortic aneurysm ----
  message("Add aortic aneurysm")

  id <- "ebi-a-GCST90018783"

  gwas[nrow(gwas) + 1, ] <- c(
    map_snps = "",
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
    ncase_col = "",
    ncase = ao[ao$id == id, ]$sample_size - ao[ao$id == id, ]$ncontrol,
    ncontrol_col = "",
    ncontrol = ao[ao$id == id, ]$ncontrol,
    chr_col = "",
    pos_col = "",
    build = "GRCh37",
    adj = "",
    ukb = ((1520 + 128425) / ao[ao$id == id, ]$sample_size),
    id = id,
    unit = "Binary",
    SD = NA
  )

  # Add arterial embolism and thrombosis ----
  message("Add arterial embolism and thrombosis")

  gwas[nrow(gwas) + 1, ] <- c(
    map_snps = "",
    source = "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90044001-GCST90045000/GCST90044016/GCST90044016_buildGRCh37.tsv.gz",
    data = "raw/GCST90044016_buildGRCh37.tsv.gz",
    phenotype = "Arterial embolism and thrombosis",
    phenotype_short = "AET",
    category = "binary",
    snp_col = "variant_id",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    eaf_col = "effect_allele_frequency",
    beta = "beta",
    se_col = "standard_error",
    pval_col = "p_value",
    samplesize_col = "N",
    samplesize = 456348,
    ncase_col = "",
    ncase = 454,
    ncontrol_col = "",
    ncontrol = 455894,
    chr_col = "chromosome",
    pos_col = "base_pair_location",
    build = "GRCh37",
    adj = "",
    ukb = 1,
    id = "",
    unit = "Binary",
    SD = NA
  )

  # Add pulmonary embolism ----
  message("Add pulmonary embolism")

  id <- "finn-b-I9_PULMEMB"

  gwas[nrow(gwas) + 1, ] <- c(
    map_snps = "",
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
    ncase_col = "",
    ncase = 4185,
    ncontrol_col = "",
    ncontrol = 214228,
    chr_col = "",
    pos_col = "",
    build = "GRCh37",
    adj = "",
    ukb = 0,
    id = id,
    unit = "Binary",
    SD = NA
  )

  # Format ----

  gwas$samplesize <- as.numeric(gwas$samplesize)
  gwas$ukb <- as.numeric(gwas$ukb)

  # Save ----

  data.table::fwrite(gwas, "data/gwas.csv")

  # Return ----

  return(gwas)
}
