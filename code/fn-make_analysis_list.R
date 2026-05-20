make_analysis_list <- function(
  gwas,
  target_traits = c("eGFRr", "eGFRcrea", "eGFRcys")
) {
  # List target traits against all available outcomes
  exp_df <- tidyr::crossing(
    exposure = target_traits,
    outcome = gwas$phenotype_short
  )

  # List target traits against all available exposures
  out_df <- tidyr::crossing(
    exposure = gwas$phenotype_short,
    outcome = target_traits
  )

  # Bind analyses lists
  df <- rbind(exp_df, out_df)

  # Remove analyses where the exposure and outcome are the same
  df <- df[df$exposure != df$outcome, ]

  # Return analyses list
  return(df)
}
