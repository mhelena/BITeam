#!/usr/bin/env Rscript
# Collects the 32 PRSice .summary files (best threshold per run) into one table.
# Usage: Rscript scripts/02_M1_summarise_PRSice.R config.sh

args <- commandArgs(trailingOnly = FALSE)
source(file.path(dirname(sub("^--file=", "", grep("^--file=", args, value = TRUE))), "config.R"))
cfg <- read_config()

SCO <- cfg("PRS_DIR")
dir.create(cfg("OUT_DIR"), showWarnings = FALSE, recursive = TRUE)
OUT <- file.path(cfg("OUT_DIR"), "PRS_summary_table_noPostSurg.txt")

gwas <- c(epilepsy = "Epilepsy", intelligence = "Intelligence",
          longevity = "Longevity", AD = "AD",
          brainage_gm = "BrainAge GM", brainage_wm = "BrainAge WM",
          brainage_gwm = "BrainAge GWM", frailty = "Frailty")
comp <- c(syndrome = "Focal vs Generalised", DRE = "DRE vs non-DRE",
          convulsions = "Convulsions vs No", currentVPA = "On VPA vs Not")

res <- list()
for (g in names(gwas)) for (c in names(comp)) {
  f <- file.path(SCO, sprintf("PRS_%s_%s.summary", g, c))
  if (!file.exists(f)) { cat("missing:", f, "\n"); next }
  s <- read.delim(f, stringsAsFactors = FALSE)
  res[[length(res)+1]] <- data.frame(
    GWAS = gwas[g], Comparison = comp[c],
    P_value = s$P, PRS_R2 = s$PRS.R2, Full_R2 = s$Full.R2, Null_R2 = s$Null.R2,
    Coefficient = s$Coefficient, SE = s$Standard.Error,
    Threshold = s$Threshold, Num_SNP = s$Num_SNP
  )
}
res <- do.call(rbind, res)
res <- res[order(res$P_value), ]
write.table(res, OUT, sep = "\t", row.names = FALSE, quote = FALSE)
cat("Saved:", OUT, "\n", nrow(res), "rows |",
    sum(res$P_value < 0.05), "raw p < 0.05\n")
