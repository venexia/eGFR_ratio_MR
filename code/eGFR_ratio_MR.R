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
message("Make analysis list")

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

# Analysis ----
message("Analysis")

# Extract instruments ----
for (i in unique(analyses$exposure)) {
  message(paste0("Make instrument for ", i))
  suppressMessages(prepare_gwas(
    gwas = gwas[gwas$phenotype_short == i, ],
    type = "exposure",
    p_threshold = 5e-8,
    clump = TRUE,
    clump_kb = 10000,
    clump_r2 = 0.001
  ))
}

# Perform MR ----
for (i in 1:nrow(analyses)) {
  message(paste0(
    "Perform MR of ",
    analyses$exposure[i],
    " on ",
    analyses$outcome[i]
  ))
  tmp <- perform_mr(
    gwas = gwas,
    exp_name = analyses$exposure[i],
    out_name = analyses$outcome[i],
    sf = TRUE
  )
}

# Combine output ----

results <- format_mr_results(filepath = "data/mr/", metadata = gwas)
data.table::fwrite(results, "output/results.csv", row.names = FALSE)
