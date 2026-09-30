# Tidy MR results ----

df_mr <- data.table::fread("output/results.csv", data.table = FALSE)

df_mr <- df_mr[
  df_mr$method %in%
    c("Inverse variance weighted", "Wald ratio") &
    df_mr$exposure == "eGFR_ratio",
  c("outcome", "est", "lci", "uci")
]

df_mr$uci <- as.numeric(df_mr$uci)

df_mr <- dplyr::rename(
  df_mr,
  "mr_trait" = "outcome",
  "mr_est" = "est",
  "mr_lci" = "lci",
  "mr_uci" = "uci"
)

# Tidy obs results ----

df_obs <- data.table::fread("raw/All_regression_results_combined_OR.csv")

df_obs$estimate_lci <- df_obs$estimate - qnorm(0.975) * df_obs$std_error
df_obs$estimate_uci <- df_obs$estimate + qnorm(0.975) * df_obs$std_error

df_obs$obs_est <- ifelse(
  df_obs$analysis_outcome_type == "Binary",
  df_obs$OR,
  ifelse(
    df_obs$analysis_outcome_type == "Continuous (Normalised)",
    df_obs$estimate,
    NA
  )
)
df_obs$obs_lci <- ifelse(
  df_obs$analysis_outcome_type == "Binary",
  df_obs$lower_CI,
  ifelse(
    df_obs$analysis_outcome_type == "Continuous (Normalised)",
    df_obs$estimate_lci,
    NA
  )
)
df_obs$obs_uci <- ifelse(
  df_obs$analysis_outcome_type == "Binary",
  df_obs$upper_CI,
  ifelse(
    df_obs$analysis_outcome_type == "Continuous (Normalised)",
    df_obs$estimate_uci,
    NA
  )
)

df_obs$obs_adj <- ifelse(
  df_obs$analysis_adjusted == "Unadjusted",
  "unadj",
  ifelse(df_obs$analysis_adjusted == "Adjusted - age and sex", "adj", "")
)

df_obs <- dplyr::rename(
  df_obs,
  "obs_trait" = "outcome_variable"
)

df_obs <- df_obs[, c("obs_trait", "obs_adj", "obs_est", "obs_lci", "obs_uci")]

df_obs <- tidyr::pivot_wider(
  df_obs,
  names_from = "obs_adj",
  values_from = c("obs_est", "obs_lci", "obs_uci")
)

# Map traits ----

map_traits <- data.table::fread("lib/label_ref.csv", data.table = FALSE)
map_traits <- dplyr::rename(map_traits, "outcome" = "trait")
outcome_id <- c("outcome", "trait_group", "type")
map_traits <- unique(map_traits[, c("obs_trait", "mr_trait", outcome_id)])

df_obs <- merge(
  df_obs,
  map_traits[, c("obs_trait", outcome_id)],
  by = "obs_trait",
  all.x = TRUE
)

df_mr <- merge(
  df_mr,
  map_traits[, c("mr_trait", outcome_id)],
  by = "mr_trait",
  all.x = TRUE
)

# Combine obs and MR results ----

df <- merge(df_obs, df_mr, by = outcome_id, all = TRUE)
df$exposure <- "eGFR ratio"

# Pivot ----

adj <- "_adj"

df <- df[, c(
  "exposure",
  outcome_id,
  paste0("obs_", c("est", "lci", "uci"), adj),
  paste0("mr_", c("est", "lci", "uci"))
)]

colnames(df) <- gsub(adj, "", colnames(df))
diff_outcome <- df[is.na(df$obs_est) | is.na(df$mr_est), ]$outcome

df <- tidyr::pivot_longer(
  df,
  cols = c(
    paste0("obs_", c("est", "lci", "uci")),
    paste0("mr_", c("est", "lci", "uci"))
  ),
  names_to = c("obs_mr", ".value"),
  names_pattern = "(obs|mr)_(est|lci|uci)"
)

df <- df[!is.na(df$est), ]
df$obs_mr <- ifelse(df$obs_mr == "obs", "Observational", df$obs_mr)
df$obs_mr <- ifelse(df$obs_mr == "mr", "MR", df$obs_mr)
df$outcome <- ifelse(
  df$outcome %in% diff_outcome,
  paste0(df$outcome, "*"),
  df$outcome
)

# Plot ----

for (i in c("binary", "continuous")) {
  p <- ggplot2::ggplot(
    df[df$type == i, ],
    mapping = ggplot2::aes(
      x = est,
      y = forcats::fct_rev(outcome),
      colour = obs_mr
    )
  )

  if (i == "binary") {
    p <- p +
      ggplot2::geom_vline(xintercept = 1, col = "dark grey") +
      ggplot2::scale_x_continuous(
        transform = "log",
        limits = c(2^-10, 2^6),
        breaks = 2^seq(-100, 100, 2),
        labels = signif(2^seq(-100, 100, 2), 2),
        name = "Beta and 95% confidence interval\nper unit increase in eGFR ratio"
      )
  }

  if (i == "continuous") {
    p <- p +
      ggplot2::geom_vline(xintercept = 0, col = "dark grey") +
      ggplot2::scale_x_continuous(
        limits = c(-2.5, 2),
        breaks = seq(-100, 100, 0.5),
        name = "Odds ratio and 95% confidence interval\nper unit increase in eGFR ratio"
      )
  }
  p +
    ggplot2::geom_linerange(
      ggplot2::aes(xmin = lci, xmax = uci),
      alpha = 0.4,
      size = 1,
      position = ggplot2::position_dodge(width = 0.5)
    ) +
    ggplot2::geom_point(
      shape = 15,
      size = 0.5,
      position = ggplot2::position_dodge(width = 0.5)
    ) +
    ggplot2::labs(
      y = "",
      colour = "Analysis"
    ) +
    ggplot2::scale_color_manual(
      values = c("#d95f02", "#7570b3"),
      breaks = c("Observational", "MR")
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      legend.position = "bottom",
      text = ggplot2::element_text(size = 8)
    ) +
    ggplot2::facet_wrap(
      ggplot2::vars(trait_group),
      scales = "free_y",
      space = "free_y"
    ) +
    ggplot2::theme(strip.text = ggplot2::element_text(angle = 0, hjust = 0))

  # Save plot
  ggplot2::ggsave(
    filename = paste0("output/compare_obs_", i, ".jpeg"),
    dpi = 300,
    width = 210,
    height = 297,
    unit = "mm",
    scale = 0.7
  )
}
