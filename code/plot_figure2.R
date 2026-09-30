library(ggplot2)
library(dplyr)
library(cowplot)

# Load observational results ----

df <- data.table::fread(
  "raw/All_regression_results_combined_OR.csv",
  data.table = FALSE
)

# Calculate CI ----

df$lci <- df$estimate - (qnorm(0.975) * df$std_error)
df$uci <- df$estimate + (qnorm(0.975) * df$std_error)

# Restrict to relevant results ----

df <- df[,
  c(
    "analysis_exposure",
    "outcome_variable",
    "analysis_adjusted",
    "estimate",
    "lci",
    "uci"
  )
]

# Rename analaysis_adjusted ----

df$analysis_adjusted <- ifelse(
  df$analysis_adjusted == "Adjusted - age and sex",
  "Age and sex adjusted",
  df$analysis_adjusted
)

df$analysis_adjusted <- factor(
  df$analysis_adjusted,
  levels = c("Unadjusted", "Age and sex adjusted")
)

# Map traits ----

map_traits <- data.table::fread("lib/label_ref.csv", data.table = FALSE)
map_traits <- dplyr::rename(map_traits, "outcome_variable" = "obs_trait")

df <- merge(
  df,
  map_traits[, c("outcome_variable", "trait", "trait_group", "type")],
  by = "outcome_variable",
  all.x = TRUE
)

# Plot

pd <- position_dodge(width = 0.5)

ggplot(
  df,
  aes(
    x = estimate,
    y = forcats::fct_rev(trait),
    shape = analysis_adjusted
  )
) +
  geom_vline(xintercept = 0, col = "grey50") +
  geom_point(position = pd, na.rm = TRUE) +
  geom_errorbarh(
    aes(xmin = lci, xmax = uci),
    position = pd,
    width = 0,
    na.rm = TRUE
  ) +
  facet_grid(
    trait_group ~ .,
    scales = "free_y",
    space = "free_y"
  ) +
  scale_x_continuous(
    limits = log(c(0.5, 1.05)),
    breaks = log(seq(0.5, 1.05, 0.1)),
    labels = round(log(seq(0.5, 1.05, 0.1)), 2),
    name = "Beta and 95% confidence interval per standard deviation",
    sec.axis = ggplot2::sec_axis(
      trans = ~ exp(.),
      name = "Odds ratio and 95% confidence interval per standard deviation",
      breaks = seq(0.5, 1.05, 0.1)
    )
  ) +
  labs(
    y = NULL,
    shape = NULL
  ) +
  theme_bw() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "bottom",
    plot.margin = margin(t = 20, r = 10, b = 10, l = 10)
  )


# Save plot
ggplot2::ggsave(
  filename = "output/figure2.jpeg",
  dpi = 300,
  width = 297,
  height = 210,
  unit = "mm",
  scale = 1.1
)
