#!/usr/bin/env Rscript
# Writes the four binary phenotype files used as PRSice contrasts (syndrome, DRE, convulsions, VPA).
# Usage: Rscript scripts/00_create_pheno_files.R config.sh

args <- commandArgs(trailingOnly = FALSE)
source(file.path(dirname(sub("^--file=", "", grep("^--file=", args, value = TRUE))), "config.R"))
cfg <- read_config()

META <- cfg("METADATA")
FAM  <- paste0(cfg("TARGET"), ".fam")
OUT  <- cfg("PHENO_DIR")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

meta <- read.delim(META, stringsAsFactors = FALSE, check.names = FALSE)
meta <- meta[!(meta$MRI_findings %in% c("Post surgery", "Large lesion")), ]
cat("Samples kept:", nrow(meta), "\n")

fam <- read.table(FAM, header = FALSE, stringsAsFactors = FALSE)
colnames(fam)[1:2] <- c("real_FID", "real_IID")

# case = 1, control = 0; any other value (missing, not known) is dropped
write_pheno <- function(col, case_val, ctrl_val, filename) {
  v <- meta[[col]]
  keep <- v %in% c(case_val, ctrl_val)
  sub  <- meta[keep, ]
  idx  <- match(sub$FID, fam$real_FID)
  out  <- data.frame(
    FID = fam$real_FID[idx],
    IID = fam$real_IID[idx],
    Pheno = ifelse(v[keep] == case_val, 1, 0)
  )
  out <- out[!is.na(out$FID), ]
  write.table(out, file.path(OUT, filename), sep = "\t",
              row.names = FALSE, quote = FALSE)
  cat(sprintf("  %-25s n=%d  (cases=%d ctrls=%d)\n",
              filename, nrow(out), sum(out$Pheno == 1), sum(out$Pheno == 0)))
}

write_pheno("EpilepsySyndrome_forPRS", "Focal", "Generalised",
            "pheno_syndrome.txt")
write_pheno("DRE", "DRE at time of scan", "NO DRE at time of scan",
            "pheno_DRE.txt")
write_pheno("HistoryOfConvulsions", "History of convulsions",
            "NO History of convulsions", "pheno_convulsions.txt")
write_pheno("VPA", "VPA exposure", "NO VPA exposure",
            "pheno_currentVPA.txt")
