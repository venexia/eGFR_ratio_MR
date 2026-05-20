mb <- 0


hpse <- readRDS("data/hpse.rds")

exp <- hpse
exp <- exp[
  exp$CHROM == 4 & exp$GENPOS >= 83292461 - mb & exp$GENPOS <= 83335153 + mb,
]

uacr <- vroom::vroom(
  "raw/formatted_20180517-UACR_overall-EA-nstud_18-SumMac_400.tbl.rsid.gz"
)

out <- merge(uacr, )
