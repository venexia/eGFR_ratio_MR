make_analysis_list <- function(
  gwas,
  target_traits = c("eGFRr", "eGFRcrea", "eGFRcys")
) {
  # All available traits
  all_traits <- unique(gwas$phenotype_short)

  # Remove target traits from the non-target pool
  non_target_traits <- setdiff(all_traits, target_traits)

  # Target traits as exposures against non-target outcomes
  exp_df <- tidyr::crossing(
    exposure = target_traits,
    outcome = non_target_traits
  )

  # Non-target exposures against target traits as outcomes
  out_df <- tidyr::crossing(
    exposure = non_target_traits,
    outcome = target_traits
  )

  # Combine analyses
  df <- dplyr::bind_rows(exp_df, out_df)

  # Safety check: remove any self-comparisons
  df <- df[df$exposure != df$outcome, ]

  return(df)
}
