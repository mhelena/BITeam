#!/usr/bin/env Rscript
# Manual review sheet: one row per carrier x variant with canonical annotation, age and HPO terms.
# Usage: source $PIPELINE_CONFIG; Rscript 07_review_sheet.R
cfg <- function(x) {
  v <- Sys.getenv(x)
  if (!nzchar(v)) stop(x, " is not set; source config.sh first", call. = FALSE)
  v
}

CARRIER_VAR <- cfg("CARRIER_VARIANTS")
VAR_ANNOT   <- cfg("VAR_ANNOT")
HPO_WIDE    <- cfg("HPO_WIDE")
OUT_FILE    <- file.path(cfg("OUT_DIR"),
                         sprintf("review_sheet_%s_%s.tsv", cfg("GENE"), format(Sys.Date())))

stopifnot(file.exists(CARRIER_VAR), file.exists(VAR_ANNOT), file.exists(HPO_WIDE))

read_tsv <- function(f) {
  x <- read.table(f, header = TRUE, sep = "\t", quote = "",
                  comment.char = "", stringsAsFactors = FALSE,
                  check.names = FALSE, na.strings = c("", "NA", "."))
  cat(sprintf("  %-45s  %6d rows  %2d cols\n", basename(f), nrow(x), ncol(x)))
  x
}

cat("Reading inputs:\n")
cv  <- read_tsv(CARRIER_VAR)
va  <- read_tsv(VAR_ANNOT)
hp  <- read_tsv(HPO_WIDE)

# Older carrier tables were written without a header.
if (!all(c("IID","CHROM","POS","REF","ALT") %in% names(cv))) {
  cat("  (no expected header in carrier table; re-reading with positional names)\n")
  cv <- read.table(CARRIER_VAR, header = FALSE, sep = "\t", quote = "",
                   comment.char = "", stringsAsFactors = FALSE)
  stopifnot(ncol(cv) >= 5)
  names(cv)[1:5] <- c("IID","CHROM","POS","REF","ALT")
}

cv$POS <- as.integer(cv$POS)
va$POS <- as.integer(va$POS)

key <- function(d) paste(d$CHROM, d$POS, d$REF, d$ALT, sep = ":")
cv$vkey <- key(cv)
va$vkey <- key(va)

va_cols <- c("vkey","Consequence","HGVSc","HGVSp","REVEL",
             "SIFT_pred","SIFT_score",
             "Polyphen2_HVAR_pred","Polyphen2_HVAR_score")
va_cols <- intersect(va_cols, names(va))
merged  <- merge(cv[, c("IID","CHROM","POS","REF","ALT","vkey")],
                 va[, va_cols],
                 by = "vkey", all.x = TRUE)

# Variants dropped by the canonical/NMD filter are flagged rather than lost; if they
# should go, add the carriers to EXCLUDE_IDS and rerun from step 4.
missing_annot <- is.na(merged$HGVSp) | merged$HGVSp == ""
if (any(missing_annot)) {
  merged$HGVSp[missing_annot] <- "<no_canonical_missense>"
  cat(sprintf("\nWARNING: %d rows have no canonical annotation (likely NMD or intronic):\n",
              sum(missing_annot)))
  print(unique(merged[missing_annot, c("CHROM","POS","REF","ALT","IID")]))
}

hp_cols <- c("IID","participant_id","year_of_birth","age_at_cut",
             "n_hpo_terms","hpo_terms")
hp_cols <- intersect(hp_cols, names(hp))
final <- merge(merged, hp[, hp_cols], by = "IID", all.x = TRUE)

# Carriers missing from the HPO table were removed by the filters or the age cut.
n_before <- nrow(final)
final <- final[!is.na(final$participant_id), ]
cat(sprintf("\nRows dropped (not in filtered, age-checked HPO table): %d\n",
            n_before - nrow(final)))

front <- c("IID","participant_id","year_of_birth","age_at_cut",
           "CHROM","POS","REF","ALT",
           "HGVSc","HGVSp","Consequence",
           "REVEL","SIFT_pred","SIFT_score",
           "Polyphen2_HVAR_pred","Polyphen2_HVAR_score",
           "n_hpo_terms","hpo_terms")
front <- intersect(front, names(final))
final <- final[, c(front, setdiff(names(final), c(front, "vkey")))]
final <- final[order(final$CHROM, final$POS, final$REF, final$ALT, final$IID), ]

write.table(final, OUT_FILE, sep = "\t", row.names = FALSE,
            quote = FALSE, na = "")
cat(sprintf("\nWrote %d rows to %s\n", nrow(final), OUT_FILE))

cat(sprintf("\nUnique samples in output:  %d\n", length(unique(final$IID))))
cat(sprintf("Unique variants in output: %d\n",
            length(unique(paste(final$CHROM, final$POS, final$REF, final$ALT)))))
cat(sprintf("Rows without HPO recorded: %d\n",
            sum(is.na(final$hpo_terms) | final$hpo_terms == "")))
if ("REVEL" %in% names(final)) {
  n_revel <- sum(!is.na(as.numeric(final$REVEL)))
  cat(sprintf("Rows with REVEL score:     %d / %d\n", n_revel, nrow(final)))
}
