library(ggplot2)
library(dplyr)
library(cowplot)

# Load MR results ----

df <- data.table::fread("output/results.csv", data.table = FALSE)

# Exclude phenotypes ----

excl_phenotypes <- c(
  "Arterial embolism and thrombosis", # no instruments
  "Chronic renal failure", # not interpretable
  "Chronic kidney disease" # not interpretable
)

df <- df[!(df$exposure %in% excl_phenotypes), ]
df <- df[!(df$outcome %in% excl_phenotypes), ]

# Restrict to relevant MR results ----

df <- df[
  df$method %in%
    c("Inverse variance weighted", "Wald ratio") &
    (df$exposure == "eGFR_ratio" | df$outcome == "eGFR_ratio"),
  c("exposure", "outcome", "b", "b_lci", "b_uci")
]

# Flip MR results ----

df$flip_b <- -1 * df$b
df$flip_b_lci <- -1 * df$b_uci
df$flip_b_uci <- -1 * df$b_lci

# Label MR results ----

df$mr <- paste0(df$exposure, " > ", df$outcome)
df$direction <- ifelse(
  df$exposure == "eGFR_ratio",
  "eGFR ratio > Phenotypes",
  ifelse(df$outcome == "eGFR_ratio", "Phenotypes > eGFR ratio", "")
)
df$phenotype <- ifelse(
  df$exposure == "eGFR_ratio",
  df$outcome,
  ifelse(df$outcome == "eGFR_ratio", df$exposure, "")
)

# Map traits ----

map_traits <- data.table::fread("lib/label_ref.csv", data.table = FALSE)
map_traits <- dplyr::rename(map_traits, "phenotype" = "mr_trait")

df <- merge(
  df,
  map_traits[, c("phenotype", "trait_group", "type")],
  by = "phenotype",
  all.x = TRUE
)

# Plot

x_min <- -3
x_max <- 3
tick <- 1

df2 <- df |>
  mutate(
    b_lci_plot = pmax(b_lci, x_min),
    b_uci_plot = pmin(b_uci, x_max),
    flip_b_lci_plot = pmax(flip_b_lci, x_min),
    flip_b_uci_plot = pmin(flip_b_uci, x_max),
    lci_trunc = b_lci < x_min,
    uci_trunc = b_uci > x_max,
    flip_lci_trunc = flip_b_lci < x_min,
    flip_uci_trunc = flip_b_uci > x_max
  )

df2 <- df2 |>
  mutate(
    est = if_else(exposure == "eGFR_ratio", flip_b, b),
    lci = if_else(exposure == "eGFR_ratio", flip_b_lci_plot, b_lci_plot),
    uci = if_else(exposure == "eGFR_ratio", flip_b_uci_plot, b_uci_plot),
    lci_trunc_flag = if_else(
      exposure == "eGFR_ratio",
      flip_lci_trunc,
      lci_trunc
    ),
    uci_trunc_flag = if_else(
      exposure == "eGFR_ratio",
      flip_uci_trunc,
      uci_trunc
    )
  )

p <- ggplot(
  df2,
  aes(
    x = est,
    y = forcats::fct_rev(phenotype)
  )
) +
  geom_vline(xintercept = 0, col = "darkred", linetype = 2) +

  geom_point(na.rm = TRUE) +

  # main (clipped) CI bars
  geom_errorbarh(
    aes(xmin = lci, xmax = uci),
    width = 0,
    na.rm = TRUE
  ) +

  # left arrow (CI goes below xmin)
  geom_segment(
    data = df2 %>% filter(lci_trunc_flag),
    aes(x = x_min, xend = x_min - 0.15, y = phenotype, yend = phenotype),
    arrow = arrow(length = unit(0.15, "cm")),
    inherit.aes = FALSE
  ) +

  # right arrow (CI goes above xmax)
  geom_segment(
    data = df2 %>% filter(uci_trunc_flag),
    aes(x = x_max, xend = x_max + 0.15, y = phenotype, yend = phenotype),
    arrow = arrow(length = unit(0.15, "cm")),
    inherit.aes = FALSE
  ) +

  facet_grid(
    trait_group ~ direction,
    scales = "free_y",
    space = "free_y"
  ) +

  scale_x_continuous(
    limits = c(x_min - 0.15, x_max + 0.15),
    breaks = log(2^seq(-1000, 1000, tick)),
    labels = round(log(2^seq(-1000, 1000, tick)), 1),
    name = "Beta and 95% confidence interval",
    sec.axis = ggplot2::sec_axis(
      trans = ~ exp(.),
      name = "Odds ratio and 95% confidence interval",
      breaks = 2^seq(-1000, 1000, tick),
      labels = as.character(c(
        formatC(2^seq(-1000, -11, tick), format = "e", digits = 1),
        round(2^seq(-10, 1000, tick), 3)
      ))
    )
  ) +

  labs(
    y = NULL
  ) +

  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(t = 20, r = 10, b = 10, l = 10)
  )

ggdraw(p) +
  draw_label(
    "Effect of an SD decrease in eGFR ratio on cardiometabolic risk factors and outcomes",
    x = 0.34,
    y = 0.99,
    size = 10,
    fontface = "bold"
  ) +
  draw_label(
    "Effect of cardiometabolic risk factors and outcomes* on eGFR ratio", # an SD increase (continuous) or unit increase in log-odds of liability (binary)
    x = 0.76,
    y = 0.99,
    size = 10,
    fontface = "bold"
  )


# Save plot
ggplot2::ggsave(
  filename = "output/figure1.jpeg",
  dpi = 300,
  width = 297,
  height = 210,
  unit = "mm",
  scale = 1.1
)
