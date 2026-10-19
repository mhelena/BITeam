#!/usr/bin/env Rscript
# Variant table with Sample1..SampleN IDs for Airlock export, plus an ID mapping that stays in the RE.
# Usage: source $PIPELINE_CONFIG; Rscript 08_airlock_table.R
cfg <- function(x) {
  v <- Sys.getenv(x)
  if (!nzchar(v)) stop(x, " is not set; source config.sh first", call. = FALSE)
  v
}

GENE       <- cfg("GENE")
CARRIERS_F <- cfg("FINAL_SAMPLES")
CV_F       <- cfg("CARRIER_VARIANTS")
VA_F       <- cfg("VAR_ANNOT")

OUT_DIR  <- cfg("AIRLOCK_DIR")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
OUT_ANON <- file.path(OUT_DIR, sprintf("%s_controls_airlock_variants.tsv", GENE))
# Kept outside the export folder so it cannot be picked up with the Airlock file.
OUT_MAP  <- file.path(cfg("OUT_DIR"), sprintf("%s_controls_sample_id_mapping.tsv", GENE))

for (f in c(CARRIERS_F, CV_F, VA_F)) {
    if (!file.exists(f)) stop(sprintf("MISSING: %s", f))
}

carriers <- readLines(CARRIERS_F)
carriers <- unique(trimws(carriers[nzchar(carriers)]))
cat(sprintf("Final %s control cohort: %d carriers\n", GENE, length(carriers)))

read_tsv <- function(f, ...) read.table(f, header = TRUE, sep = "\t", quote = "",
                                         comment.char = "", stringsAsFactors = FALSE,
                                         check.names = FALSE, na.strings = c("","NA","."), ...)
cv <- tryCatch(read_tsv(CV_F), error = function(e) NULL)
if (is.null(cv) || !all(c("IID","CHROM","POS","REF","ALT") %in% names(cv))) {
    cv <- read.table(CV_F, header = FALSE, sep = "\t", quote = "",
                     comment.char = "", stringsAsFactors = FALSE)
    stopifnot(ncol(cv) >= 5)
    names(cv)[1:5] <- c("IID","CHROM","POS","REF","ALT")
}
cv$POS <- as.integer(cv$POS)

cv <- cv[cv$IID %in% carriers, ]
cat(sprintf("Carrier-variant rows for the final list: %d\n", nrow(cv)))

va <- read_tsv(VA_F)
va$POS <- as.integer(va$POS)

key <- function(d) paste(d$CHROM, d$POS, d$REF, d$ALT, sep = ":")
cv$vkey <- key(cv); va$vkey <- key(va)

wanted_va_cols <- intersect(
    c("vkey", "HGVSc", "HGVSp", "Consequence", "REVEL",
      "SIFT_pred", "Polyphen2_HVAR_pred", "CANONICAL", "MANE_SELECT",
      "SYMBOL", "BIOTYPE"),
    names(va))
merged <- merge(cv, va[, wanted_va_cols, drop = FALSE],
                by = "vkey", all.x = TRUE, sort = FALSE)

# Sorting by IID first makes the Sample numbering reproducible between runs.
iid_order <- sort(unique(merged$IID))
sample_map <- setNames(sprintf("Sample%d", seq_along(iid_order)), iid_order)
merged$Sample_ID <- sample_map[merged$IID]

participant_col <- intersect(c("participant_id","participant_ID","PID"), names(merged))
if (length(participant_col) == 1) {
    map_out <- unique(merged[, c("Sample_ID","IID", participant_col)])
    names(map_out)[3] <- "participant_id"
} else {
    map_out <- unique(merged[, c("Sample_ID","IID")])
    map_out$participant_id <- NA_character_
}
map_out <- map_out[order(as.integer(sub("^Sample","", map_out$Sample_ID))), ]

anon_cols <- c("Sample_ID","CHROM","POS","REF","ALT",
               "HGVSc","HGVSp","Consequence","REVEL",
               "SIFT_pred","Polyphen2_HVAR_pred","CANONICAL","MANE_SELECT")
anon_cols <- intersect(anon_cols, names(merged))
anon <- merged[, anon_cols, drop = FALSE]
anon <- anon[order(as.integer(sub("^Sample","", anon$Sample_ID)), anon$POS), ]

write.table(anon,    OUT_ANON, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
write.table(map_out, OUT_MAP,  sep = "\t", row.names = FALSE, quote = FALSE, na = "")

cat(sprintf("\nWrote %s\n", OUT_ANON))
cat(sprintf("Wrote %s   (internal, do NOT export)\n", OUT_MAP))
cat(sprintf("N samples in anonymised table : %d\n", length(unique(anon$Sample_ID))))
cat(sprintf("N variant rows                : %d\n", nrow(anon)))

cat("\nFirst 10 rows of the anonymised table:\n")
print(head(anon, 10), row.names = FALSE)

cat("\nConsequence breakdown:\n")
if ("Consequence" %in% names(anon))
    print(table(anon$Consequence), row.names = FALSE)

per_sample <- table(anon$Sample_ID)
if (any(per_sample > 1)) {
    warning("Some samples have more than one ", GENE, " variant row: either several\n",
            "  missense variants in one carrier (worth reviewing) or duplicate rows from the join.\n",
            "  Samples affected: ",
            paste(names(per_sample[per_sample > 1]), collapse = ", "))
}
