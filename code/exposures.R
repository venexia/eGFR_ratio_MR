exp <- data.frame(
  SNP = character(),
  chr.exposure = character(),
  pos.exposure = character(),
  effect_allele.exposure = character(),
  other_allele.exposure = character(),
  eaf.exposure = numeric(),
  beta.exposure = numeric(),
  se.exposure = numeric(),
  pval.exposure = numeric(),
  samplesize.exposure = numeric(),
  ncase.exposure = numeric(),
  ncontrol.exposure = numeric(),
  exposure = character(),
  mr_keep.exposure = character(),
  pval_origin.exposure = character()
)

for (i in 1:nrow(gwas)) {
  message(paste0(
    "Exposure ",
    i,
    " of ",
    nrow(gwas),
    ": ",
    gwas[i, "phenotype"]
  ))

  if (
    file.exists(paste0("data/exposure_", gwas[i, "phenotype_short"], ".csv"))
  ) {
    message("Skipped. File already exists.")
  } else {
    tmp <- suppressMessages(prepare_gwas(
      map_snps = gwas[i, "map_snps"],
      type = "exposure",
      data = gwas[i, "data"],
      phenotype = gwas[i, "phenotype"],
      phenotype_short = gwas[i, "phenotype_short"],
      category = gwas[i, "category"],
      snp_col = gwas[i, "snp_col"],
      effect_allele_col = gwas[i, "effect_allele_col"],
      other_allele_col = gwas[i, "other_allele_col"],
      eaf_col = gwas[i, "eaf_col"],
      beta_col = gwas[i, "beta_col"],
      se_col = gwas[i, "se_col"],
      pval_col = gwas[i, "pval_col"],
      samplesize_col = gwas[i, "samplesize_col"],
      samplesize = gwas[i, "samplesize"],
      ncase_col = gwas[i, "ncase_col"],
      ncase = gwas[i, "ncase"],
      ncontrol_col = gwas[i, "ncontrol_col"],
      ncontrol = gwas[i, "ncontrol"],
      chr_col = gwas[i, "chr_col"],
      pos_col = gwas[i, "pos_col"],
      id = gwas[i, "id"],
      p_threshold = 5e-8,
      clump = TRUE,
      clump_kb = 10000,
      clump_r2 = 0.001,
      SD = gwas[i, "SD"]
    ))

    message(paste0(
      "Saving exposure data/exposure_",
      gwas[i, "phenotype_short"],
      ".csv"
    ))

    if (!is.null(nrow(tmp))) {
      exp <- plyr::rbind.fill(exp, tmp)

      data.table::fwrite(
        exp,
        paste0("data/exposure_", gwas[i, "phenotype_short"], ".csv"),
        row.names = FALSE
      )
    }
  }
}
