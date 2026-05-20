make_plot <- function(
  df,
  pheno_col = "outcome",
  est_col = "est",
  lci_col = "lci",
  uci_col = "uci",
  bin_cont = "both",
  colour = "exposure",
  plot_min,
  plot_max,
  tick = 1,
  plot_name = "plot",
  portrait = TRUE
) {
  # Rename columns ----
  df <- dplyr::rename(
    df,
    "pheno" = pheno_col,
    "est" = est_col,
    "lci" = lci_col,
    "uci" = uci_col
  )

  # Clamp CIs ----
  df <- df %>%
    dplyr::mutate(
      beyond_lci = lci < plot_min,
      beyond_uci = uci > plot_max,
      lci_clamped = pmax(lci, plot_min),
      uci_clamped = pmin(uci, plot_max)
    )

  # Base plot ----
  p <- ggplot2::ggplot(
    df,
    mapping = ggplot2::aes(
      x = est,
      y = forcats::fct_rev(pheno),
      colour = .data[[colour]]
    )
  ) +
    ggplot2::geom_vline(xintercept = 0, col = "dark grey") +
    ggplot2::geom_linerange(
      ggplot2::aes(xmin = lci_clamped, xmax = uci_clamped),
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
      y = ""
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      legend.position = "bottom",
      text = ggplot2::element_text(size = 8)
    )

  # Add x-axis scale depending on j
  if (bin_cont == "cont") {
    p <- p +
      ggplot2::scale_x_continuous(
        limits = c(plot_min, plot_max),
        breaks = seq(-10000, 10000, tick),
        name = "Beta and 95% confidence interval",
      )
  } else if (bin_cont == "bin") {
    p <- p +
      ggplot2::scale_x_continuous(
        limits = c(
          exp(plot_min),
          exp(plot_max)
        ),
        name = "Odds ratio and 95% confidence interval per one unit decrease",
        breaks = 2^seq(-1000, 1000, tick),
        labels = as.character(c(
          formatC(2^seq(-1000, -11, tick), format = "e", digits = 1),
          round(2^seq(-10, 1000, tick), 3)
        ))
      )
  } else {
    p <- p +
      ggplot2::scale_x_continuous(
        limits = c(plot_min, plot_max),
        breaks = log(2^seq(-1000, 1000, tick)),
        labels = round(log(2^seq(-1000, 1000, tick)), 1),
        name = "Beta and 95% confidence interval per one unit decrease",
        sec.axis = ggplot2::sec_axis(
          trans = ~ exp(.),
          name = "Odds ratio and 95% confidence interval per one unit decrease",
          breaks = 2^seq(-1000, 1000, tick),
          labels = as.character(c(
            formatC(2^seq(-1000, -11, tick), format = "e", digits = 1),
            round(2^seq(-10, 1000, tick), 3)
          ))
        )
      )
  }

  # Save plot
  ggplot2::ggsave(
    filename = paste0("output/", plot_name, ".jpeg"),
    plot = p,
    dpi = 300,
    width = ifelse(isTRUE(portrait), 210, 297),
    height = ifelse(isTRUE(portrait), 297, 210),
    unit = "mm",
    scale = 0.7
  )

  return(p)
}
