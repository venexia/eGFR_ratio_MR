batch_extract_outcome <- function(snp_list, batch_size = 50, outcomes) {
  # Batch the snp_list
  snp_batches <- split(snp_list, ceiling(seq_along(snp_list) / batch_size))

  # Loop through batches
  outcome_data_list <- lapply(snp_batches, function(snps) {
    TwoSampleMR::extract_outcome_data(snps, outcomes = outcomes)
  })

  # Combine all batches
  outcome_data <- dplyr::bind_rows(outcome_data_list)

  # Return outcome data
  return(outcome_data)
}
