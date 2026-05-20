# List instrument SNPs ----
message("List instrument SNPs")

snplist <- data.frame(snp = character())

for (i in list.files("data", pattern = "exposure_", full.names = TRUE)) {
  tmp <- data.table::fread(i, select = c("SNP"), data.table = FALSE)
  snplist <- rbind(snplist, tmp)
}

snplist <- unique(snplist)
data.table::fwrite(snplist, "data/snplist.csv", row.names = FALSE)

out <- data.frame(
  SNP = character(),
  chr.outcome = character(),
  pos.outcome = character(),
  effect_allele.outcome = character(),
  other_allele.outcome = character(),
  eaf.outcome = numeric(),
  beta.outcome = numeric(),
  se.outcome = numeric(),
  pval.outcome = numeric(),
  samplesize.outcome = numeric(),
  ncase.outcome = numeric(),
  ncontrol.outcome = numeric(),
  outcome = character()
)

for (i in 1) {
  message(paste0(
    "Preparing outcome ",
    i,
    " of ",
    nrow(gwas),
    ": ",
    gwas[i, "phenotype"]
  ))

  if (
    file.exists(paste0("data/outcome_", gwas[i, "phenotype_short"], ".csv"))
  ) {
    message("Skipped. File already exists.")
  } else {
    tmp <- suppressMessages(prepare_gwas(
      map_snps = gwas[i, "map_snps"],
      type = "outcome",
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
      instrument = unique(snplist$SNP),
      SD = gwas[i, "SD"]
    ))

    message(paste0(
      "Saving exposure data/outcome_",
      gwas[i, "phenotype_short"],
      ".csv"
    ))

    if (!is.null(nrow(tmp))) {
      out <- plyr::rbind.fill(out, tmp)

      data.table::fwrite(
        out,
        paste0("data/outcome_", gwas[i, "phenotype_short"], ".csv"),
        row.names = FALSE
      )
    }
  }
}
