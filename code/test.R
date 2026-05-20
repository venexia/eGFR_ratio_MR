source("code/fn-make_plot.R")

df <- data.table::fread("output/results.csv")

map_outcomes <- data.table::fread("lib/label_ref.csv")

df <- merge(
  df,
  map_outcomes,
  by.x = "outcome",
  by.y = "mr_trait",
  all.x = TRUE
)

df <- dplyr::rename(
  df,
  "outcome_full" = "trait",
  "outcome_group" = "trait_group"
)

egfr_exp <- df[
  df$exposure %in%
    c("eGFR_ratio", "eGFR (creatinine)", "eGFR (cystatin C)") &
    !(df$outcome %in%
      c("eGFR_ratio", "eGFR (creatinine)", "eGFR (cystatin C)")) &
    !(df$outcome_group %in% c("Renal")) &
    df$method %in% c("Inverse variance weighted", "Wald ratio"),
  c(
    "exposure",
    "outcome_full",
    "outcome_group",
    "category",
    "b",
    "b_lci",
    "b_uci"
  )
]

make_plot(
  df = egfr_exp,
  pheno_col = "outcome_full",
  est_col = "b",
  lci_col = "b_lci",
  uci_col = "b_uci",
  tick = 8,
  plot_min = -60,
  plot_max = 30,
  plot_name = "test"
)
