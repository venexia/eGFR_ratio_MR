perform_mr <- function(gwas = gwas, exp_name, out_name, sf = TRUE) {
  # Check if file already exists
  message('Check if file already exists')
  mr_filename <- paste0("data/mr/", exp_name, "_", out_name, ".csv")

  if (file.exists(mr_filename)) {
    message("Skipped. File already exists.")
    mr <- data.table::fread(mr_filename, data.table = FALSE)
  } else {
    # Get exposure data ----
    message("Get exposure data")

    exp <- data.table::fread(
      paste0("data/exposure/", exp_name, ".csv"),
      data.table = FALSE
    )

    if (nrow(exp) == 0) {
      mr <- data.frame(message = "Instrument file contains 0 rows.")
      message(mr$message)
    } else {
      # Get outcome data ----
      message("Get outcome data")

      out <- prepare_gwas(
        gwas = gwas[gwas$phenotype_short == out_name, ],
        type = "outcome",
        instrument = exp$SNP,
        save = FALSE
      )

      if (nrow(out) > 0) {
        # Harmonize data ----
        message("Harmonize data")

        dat <- TwoSampleMR::harmonise_data(
          exposure_dat = exp,
          outcome_dat = out
        )

        # Perform MR ----
        message("Perform MR")

        mr <- TwoSampleMR::mr(dat)

        # Optional: Steiger filtered MR ----

        if (isTRUE(sf)) {
          # Perform Steiger filtering ----
          message("Perform Steiger filtering")

          dat_sf <- dat
          dat_sf$sd.exposure <- 0.2
          dat_sf <- TwoSampleMR::steiger_filtering(dat)

          if (nrow(dat_sf) == 0) {
            message("All SNPs removed by Steiger filtering")
            mr_sf <- unique(mr[, c(
              "id.exposure",
              "id.outcome",
              "exposure",
              "outcome"
            )])
            mr_sf <- transform(
              mr_sf,
              method = "",
              nsnp = 0,
              b = NA_real_,
              se = NA_real_,
              pval = NA_real_
            )
          } else if (nrow(dat_sf) == nrow(dat)) {
            message("Nothing removed by Steiger filtering")
            mr_sf <- mr[
              mr$method %in% c("Inverse variance weighted", "Wald ratio"),
            ]
          } else {
            message("Run MR with Steiger filtering")
            mr_sf <- TwoSampleMR::mr(
              dat_sf,
              method_list = c("mr_wald", "mr_ivw")
            )
          }

          # Add Steiger filtering results to main results ----
          message("Add Steiger filtering results to main results")

          mr_sf$method <- paste0(mr_sf$method, "_steiger")
          mr <- rbind(mr, mr_sf)
        }

        # Record additional information ----
        message("Record additional information")

        mr$nsnp.exposure <- nrow(exp)
        mr$nsnp.outcome <- nrow(out)
      } else {
        mr <- data.frame(message = "No instrument SNPs in outcome data.")
        message(mr$message)
      }
    }

    # Save results ----
    message("Save results")

    data.table::fwrite(
      mr,
      mr_filename,
      row.names = FALSE
    )
  }

  # Return results ----
  message("Return results")

  return(mr)
}
