# Load data ----

df <- vroom::vroom("output/results.csv")
gwas <- vroom::vroom("data/gwas.csv")

#

df$outcome <- gsub(" \\|.*", "", df$outcome)
df$outcome <- gsub("body mass", "Body mass", df$outcome)
df$outcome <- gsub("HbA1C", "HbA1c", df$outcome)
df <- merge(
  df,
  gwas[, c("phenotype", "category")],
  by.x = "outcome",
  by.y = "phenotype",
  all.x = TRUE
)

df$est <- ifelse(df$category == "binary", df$or, log(df$or))
df$est_lci <- ifelse(df$category == "binary", df$lci, log(df$lci))
df$est_uci <- ifelse(df$category == "binary", df$uci, log(df$uci))

# Plot ---

ggplot2::ggplot(
  df[df$method == "Inverse variance weighted" & df$category == "continuous", ],
  mapping = ggplot2::aes(x = est, y = outcome, colour = exposure)
) +
  ggplot2::geom_vline(xintercept = 0, col = "dark grey") +
  ggplot2::geom_linerange(
    ggplot2::aes(xmin = est_lci, xmax = est_uci),
    alpha = 0.5,
    size = 1,
    position = ggplot2::position_dodge(width = 0.5)
  ) +
  ggplot2::geom_point(
    shape = 15,
    size = 0.5,
    position = ggplot2::position_dodge(width = 0.5)
  ) +
  ggplot2::labs(x = "Estimate and 95% confidence interval", y = "") +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_blank(),
    text = ggplot2::element_text(size = 8)
  )

ggplot2::ggsave(
  filename = paste0("output/outcomes_continuous.jpeg"),
  dpi = 300,
  width = 210,
  height = 297,
  unit = "mm",
  scale = 0.6
)

ggplot2::ggplot(
  df[df$method == "Inverse variance weighted" & df$category == "binary", ],
  mapping = ggplot2::aes(x = est, y = outcome, colour = exposure)
) +
  ggplot2::geom_vline(xintercept = 1, col = "dark grey") +
  ggplot2::geom_linerange(
    ggplot2::aes(xmin = est_lci, xmax = est_uci),
    alpha = 0.5,
    size = 1,
    position = ggplot2::position_dodge(width = 0.5)
  ) +
  ggplot2::geom_point(
    shape = 15,
    size = 0.5,
    position = ggplot2::position_dodge(width = 0.5)
  ) +
  ggplot2::scale_x_continuous(
    trans = "log",
    lim = c(2^-17, 2^7),
    breaks = 2^seq(-20, 8, 4)
  ) +
  ggplot2::labs(x = "Estimate and 95% confidence interval", y = "") +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_blank(),
    text = ggplot2::element_text(size = 8)
  )

ggplot2::ggsave(
  filename = paste0("output/outcomes_binary.jpeg"),
  dpi = 300,
  width = 297,
  height = 210,
  unit = "mm",
  scale = 0.6
)
