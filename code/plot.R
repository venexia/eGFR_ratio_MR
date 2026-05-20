library(magrittr)

source("code/utility.R")

# Load data ----

df <- vroom::vroom("output/results.csv")

# Flip results to represent a decrease in eGFR ratio ----

df$flip_est <- ifelse(df$category == "continuous", -1 * df$b, 1 / exp(df$b))

df$flip_lci <- ifelse(
  df$category == "continuous",
  -1 * df$b_lci,
  1 / exp(df$b_lci)
)

df$flip_uci <- ifelse(
  df$category == "continuous",
  -1 * df$b_uci,
  1 / exp(df$b_uci)
)

# Define focus ----

label_ref <- data.table::fread("lib/label_ref.csv", data.table = FALSE)
focus <- c("eGFR_ratio", "eGFR (cystatin C)", "eGFR (creatinine)")
method <- c("Inverse variance weighted")

for (i in focus) {
  # Make dataset where focus is exposure
  tmp_exp <- df[
    df$method %in%
      method &
      df$exposure == i &
      !(df$outcome %in% focus),
  ]
  if (nrow(tmp_exp) > 0) {
    tmp_exp$focus <- i
    tmp_exp$mr_trait <- tmp_exp$outcome
    tmp_exp$exp_out <- "exposure"
  }
  # Make dataset where focus is outcome
  tmp_out <- df[
    df$method %in%
      method &
      df$outcome == i &
      !(df$exposure %in% focus),
  ]
  if (nrow(tmp_out) > 0) {
    tmp_out$focus <- i
    tmp_out$mr_trait <- tmp_out$exposure
    tmp_out$exp_out <- "outcome"
  }
  tmp <- rbind(tmp_exp, tmp_out)

  # Map traits ----

  map_traits <- label_ref[
    label_ref$mr_trait != "",
    c("mr_trait", "trait", "trait_group")
  ]

  tmp <- merge(
    tmp,
    map_traits,
    by = "mr_trait",
    all.x = TRUE
  )

  tmp$phenotype <- tmp$trait

  assign(i, tmp)
}

# Define common plot elements ----

for (i in focus) {
  df_plot <- get(i)

  # Sort labels ----

  df_plot$print_est <- paste0(
    df_plot$phenotype,
    "\n",
    ifelse(df_plot$category == "binary", "OR: ", "Beta: "),
    display(df_plot$flip_est),
    " (95% CI: ",
    display(df_plot$flip_lci),
    " to ",
    display(df_plot$flip_uci),
    ")"
  )

  plot_min <- -2.5
  plot_max <- 2.5

  df_plot <- df_plot %>%
    dplyr::mutate(
      beyond_lci = b_lci < plot_min,
      beyond_uci = b_uci > plot_max,
      b_lci_clamped = pmax(b_lci, plot_min),
      b_uci_clamped = pmin(b_uci, plot_max)
    )

  for (j in unique(df_plot$exp_out)) {
    # Base plot
    p <- ggplot2::ggplot(
      df_plot[df_plot$exp_out == j, ],
      mapping = ggplot2::aes(
        x = b,
        y = forcats::fct_rev(phenotype),
        colour = exp_out
      )
    ) +
      ggplot2::geom_vline(xintercept = 0, col = "dark grey") +
      ggplot2::geom_linerange(
        ggplot2::aes(xmin = b_lci_clamped, xmax = b_uci_clamped),
        alpha = 0.5,
        size = 1,
        colour = "#7570b3",
        position = ggplot2::position_dodge(width = 0.5)
      ) +
      ggplot2::geom_point(
        shape = 15,
        size = 0.5,
        colour = "#7570b3",
        position = ggplot2::position_dodge(width = 0.5)
      ) +
      ggplot2::labs(
        y = "",
        colour = ""
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank(),
        panel.grid.major.y = ggplot2::element_blank(),
        legend.position = "bottom",
        text = ggplot2::element_text(size = 8)
      )

    # Add x-axis scale depending on j
    if (j == "exposure") {
      p <- p +
        ggplot2::scale_x_continuous(
          limits = c(plot_min, plot_max),
          breaks = log(2^seq(-100, 100)),
          labels = round(log(2^seq(-100, 100)), 1),
          name = "Beta and 95% confidence interval per one unit decrease",
          sec.axis = ggplot2::sec_axis(
            trans = ~ exp(.),
            name = "Odds ratio and 95% confidence interval per one unit decrease",
            breaks = 2^seq(-100, 100),
            labels = as.character(c(
              formatC(2^seq(-100, -11), format = "e", digits = 1),
              round(2^seq(-10, 100), 3)
            ))
          )
        )
    } else {
      p <- p +
        ggplot2::scale_x_continuous(
          limits = c(-0.35, 0.35),
          breaks = seq(-10, 10, 0.1),
          name = "Beta and 95% confidence interval",
        )
    }

    p <- p +
      ggplot2::facet_wrap(
        ggplot2::vars(trait_group),
        scales = "free_y",
        space = "free_y"
      ) +
      ggplot2::theme(strip.text = ggplot2::element_text(angle = 0, hjust = 0))

    # Save plot
    ggplot2::ggsave(
      filename = paste0("output/", j, "-", i, ".jpeg"),
      plot = p,
      dpi = 300,
      width = 210,
      height = 297,
      unit = "mm",
      scale = 0.7
    )
  }
}
