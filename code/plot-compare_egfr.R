df_full <- data.table::fread("output/results.csv")

map_outcomes <- data.table::fread("lib/label_ref.csv")

for (i in c("exposure", "outcome")) {
  df <- df_full

  if (i == "exposure") {
    df <- df[
      grepl("eGFR", df$exposure) &
        df$method == "Inverse variance weighted" &
        !grepl("eGFR", df$outcome),
    ]
    df$egfr_phenotype <- df$exposure
    df$phenotype <- df$outcome
    lim_u <- 30
    lim_l <- -110
    tick <- 10
  }

  if (i == "outcome") {
    df <- df[
      grepl("eGFR", df$outcome) &
        df$method == "Inverse variance weighted" &
        !grepl("eGFR", df$exposure),
    ]
    df$egfr_phenotype <- df$outcome
    df$phenotype <- df$exposure
    lim_u <- 0.3
    lim_l <- -0.3
    tick <- 0.1
  }

  df <- merge(
    df,
    map_outcomes,
    by.x = "phenotype",
    by.y = "mr_trait",
    all.x = TRUE
  )

  ggplot2::ggplot(
    df,
    mapping = ggplot2::aes(
      x = b,
      y = forcats::fct_rev(trait),
      colour = egfr_phenotype
    )
  ) +
    ggplot2::geom_vline(xintercept = 0, col = "dark grey") +
    ggplot2::geom_linerange(
      ggplot2::aes(xmin = b_lci, xmax = b_uci),
      alpha = 0.5,
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
      colour = ""
    ) +
    ggplot2::scale_color_manual(
      values = c("#1b9e77", "#e7298a", "#7570b3"),
      breaks = c("eGFR (cystatin C)", "eGFR (creatinine)", "eGFR_ratio")
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      legend.position = "bottom",
      text = ggplot2::element_text(size = 8)
    ) +
    ggplot2::scale_x_continuous(
      limits = c(lim_l, lim_u),
      breaks = seq(-10000, 10000, tick),
      name = ifelse(
        i == "exposure",
        "Beta and 95% confidence interval\n\n[Exposures: eGFR phenotypes (unit increase);\nOutcomes: various (SD increase)]",
        "Beta and 95% confidence interval\n\n[Exposures: various (SD increase);\nOutcomes: eGFR phenotypes (unit increase)]"
      )
    ) +
    ggplot2::facet_wrap(
      ggplot2::vars(trait_group),
      ncol = 1,
      scales = "free_y",
      space = "free_y"
    ) +
    ggplot2::theme(strip.text = ggplot2::element_text(angle = 0, hjust = 0))

  # Save plot
  ggplot2::ggsave(
    filename = paste0("output/compare_egfr_", i, ".jpeg"),
    dpi = 300,
    width = 210,
    height = 297,
    unit = "mm",
    scale = 0.7
  )
}
