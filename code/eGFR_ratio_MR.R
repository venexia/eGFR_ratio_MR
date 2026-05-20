# Source functions ----
message("Source functions")

lapply(
  list.files("code", full.names = TRUE, pattern = "fn-"),
  source
)

# Make GWAS list ----
message("Make GWAS list")

if (file.exists("data/gwas.csv")) {
  gwas <- data.table::fread("data/gwas.csv", data.table = FALSE)
} else {
  gwas <- get_gwas_info()
}

# Make analysis list ----
message("Make GWAS list")

analyses <- make_analysis_list(
  gwas = gwas,
  target_traits = c("eGFRr", "eGFRcrea", "eGFRcys")
)

# Make SNP map ----
message("Make SNP map")

if (file.exists("raw/SNPmap.csv")) {
  snp_map <- data.table::fread("raw/SNPmap.csv", data.table = FALSE)
} else {
  snp_map <- make_snp_map()
}

# Create directories ----
message("Create directories")

dir.create("data/exposure", recursive = TRUE, showWarnings = FALSE)
dir.create("data/mr", recursive = TRUE, showWarnings = FALSE)

# Prepare for analysis ----
message("Prepare for analysis")

j = unique(analyses$outcome)[1]

for (i in unique(analyses$exposure)[1:3]) {
  message(paste0("Make instrument for ", i))
  suppressMessages(prepare_gwas(
    gwas = gwas[gwas$phenotype_short == i, ],
    type = "exposure",
    p_threshold = 5e-8,
    clump = TRUE,
    clump_kb = 10000,
    clump_r2 = 0.001
  ))

  message(paste0("Perform MR of ", i, " on ", j))
  tmp <- perform_mr(
    gwas = gwas,
    exp_name = i,
    out_name = j,
    sf = TRUE
  )
}
