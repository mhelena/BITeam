#!/usr/bin/env Rscript
# Per-sample HPO terms (long and wide tables) for a sample list, after the age filter.
# Usage: source $PIPELINE_CONFIG; Rscript 05_hpo_per_sample.R [id_list]   (default: $CONTROLS_LIST)
suppressPackageStartupMessages(library(Rlabkey))

cfg <- function(x) {
  v <- Sys.getenv(x)
  if (!nzchar(v)) stop(x, " is not set; source config.sh first", call. = FALSE)
  v
}

args    <- commandArgs(trailingOnly = TRUE)
IN_FILE <- if (length(args) >= 1) args[1] else cfg("CONTROLS_LIST")
OUT_DIR <- dirname(IN_FILE)
TAG     <- tools::file_path_sans_ext(basename(IN_FILE))

# Phenotyping is thin in very young children, so carriers must be at least MIN_AGE
# at the data cut year.
MIN_AGE       <- as.integer(cfg("MIN_AGE"))
DATA_CUT_YEAR <- as.integer(cfg("DATA_CUT_YEAR"))
MAX_YOB       <- DATA_CUT_YEAR - MIN_AGE

stopifnot("input file not found" = file.exists(IN_FILE))

raw_lines <- trimws(readLines(IN_FILE))
raw_lines <- raw_lines[nchar(raw_lines) > 0]
cat(sprintf("[%s] input rows: %d\n", format(Sys.time()), length(raw_lines)))
cat(sprintf("  first line: '%s'\n", raw_lines[1]))
parts  <- strsplit(raw_lines, "\\s+")
n_cols <- max(sapply(parts, length))

fetch <- function(q) {
  cat(sprintf("  loading %s\n", q))
  labkey.selectRows(baseUrl = cfg("LABKEY_BASE"), folderPath = cfg("LABKEY_FOLDER"),
                    schemaName = "lists", queryName = q,
                    colNameOpt = "rname", maxRows = 1e7)
}
anc  <- fetch("aggregate_gvcf_sample_stats")
hpo  <- fetch("rare_diseases_participant_phenotype")
part <- fetch("participant")
yob_col <- intersect(c("year_of_birth","yearofbirth"), names(part))[1]
stopifnot("No year_of_birth column in participant table" = !is.na(yob_col))
pid_to_yob <- setNames(as.integer(part[[yob_col]]),
                       as.character(part$participant_id))

PK <- intersect(c("platekey","plate_key","sample_id"), names(anc))[1]
pk_set  <- unique(as.character(anc[[PK]]))
pid_set <- unique(as.character(anc$participant_id))

best <- list(col = NA, kind = NA, hits = 0, ids = character(0))
for (col in seq_len(n_cols)) {
  ids_col <- unique(na.omit(sapply(parts, function(x)
                                   if (length(x) >= col) x[col] else NA)))
  h_pk  <- sum(ids_col %in% pk_set)
  h_pid <- sum(ids_col %in% pid_set)
  cat(sprintf("  col %d  ->  platekey=%d  participant_id=%d\n",
              col, h_pk, h_pid))
  if (h_pk  >= max(best$hits, h_pid))
    best <- list(col=col, kind="platekey",       hits=h_pk,  ids=ids_col)
  if (h_pid >  best$hits)
    best <- list(col=col, kind="participant_id", hits=h_pid, ids=ids_col)
}
if (best$hits == 0) {
  cat("\nNo matches found. Check input format.\n"); quit(status=1)
}
cat(sprintf("\n  -> using column %d as %s  (%d matches)\n",
            best$col, best$kind, best$hits))

pk_to_pid <- setNames(as.character(anc$participant_id), as.character(anc[[PK]]))
pid_to_pk <- setNames(as.character(anc[[PK]]), as.character(anc$participant_id))

if (best$kind == "platekey") {
  matched <- best$ids[best$ids %in% pk_set]
  iid_pid <- data.frame(IID            = matched,
                        participant_id = unname(pk_to_pid[matched]),
                        stringsAsFactors = FALSE)
} else {
  matched <- best$ids[best$ids %in% pid_set]
  iid_pid <- data.frame(IID            = unname(pid_to_pk[matched]),
                        participant_id = matched,
                        stringsAsFactors = FALSE)
}
iid_pid <- iid_pid[!is.na(iid_pid$participant_id) & !is.na(iid_pid$IID), ]
iid_pid <- unique(iid_pid)
cat(sprintf("  unique samples joined: %d\n", nrow(iid_pid)))

iid_pid$year_of_birth <- pid_to_yob[iid_pid$participant_id]
iid_pid$age_at_cut    <- DATA_CUT_YEAR - iid_pid$year_of_birth
n_before <- nrow(iid_pid)
missing_yob <- sum(is.na(iid_pid$year_of_birth))
iid_pid <- iid_pid[!is.na(iid_pid$year_of_birth) &
                   iid_pid$year_of_birth <= MAX_YOB, ]
cat(sprintf("  age filter (>= %d yr at %d, YoB <= %d): %d -> %d kept  (%d dropped for missing YoB)\n",
            MIN_AGE, DATA_CUT_YEAR, MAX_YOB, n_before, nrow(iid_pid), missing_yob))

hpo$participant_id <- as.character(hpo$participant_id)
hpo_present_flag <- hpo$hpo_present %in% c(TRUE,"Yes","yes","Y",1,"1")
hpo_sub <- hpo[hpo$participant_id %in% iid_pid$participant_id, ]

hpo_cols <- c("participant_id","hpo_term","hpo_id","normalised_hpo_term",
              "hpo_present")
keep_cols <- intersect(hpo_cols, names(hpo_sub))
hpo_sub <- hpo_sub[, keep_cols]
hpo_sub$hpo_present_flag <- hpo_sub$participant_id %in%
                            hpo$participant_id[hpo_present_flag]

long <- merge(iid_pid, hpo_sub, by = "participant_id", all.x = TRUE)
long <- long[order(long$IID, long$participant_id, long$hpo_term), ]
front <- c("IID","participant_id")
long <- long[, c(front, setdiff(names(long), front))]

out_long <- file.path(OUT_DIR, sprintf("%s_hpo_long.tsv", TAG))
write.table(long, out_long, sep = "\t", row.names = FALSE,
            quote = FALSE, na = "")
cat(sprintf("\nWrote long table (%d rows): %s\n", nrow(long), out_long))

joiner <- function(v) {
  v <- v[!is.na(v) & nchar(v) > 0]
  if (!length(v)) return("")
  paste(unique(v), collapse = "; ")
}

term_col <- intersect(c("normalised_hpo_term","hpo_term"), names(long))[1]

wide <- do.call(rbind, lapply(split(long, long$IID), function(d) {
  data.frame(
    IID            = d$IID[1],
    participant_id = d$participant_id[1],
    year_of_birth  = pid_to_yob[d$participant_id[1]],
    age_at_cut     = DATA_CUT_YEAR - pid_to_yob[d$participant_id[1]],
    n_hpo_terms    = sum(!is.na(d[[term_col]]) & nchar(d[[term_col]]) > 0),
    hpo_terms      = joiner(d[[term_col]]),
    stringsAsFactors = FALSE
  )
}))
wide <- wide[order(-wide$n_hpo_terms, wide$IID), ]

out_wide <- file.path(OUT_DIR, sprintf("%s_hpo_wide.tsv", TAG))
write.table(wide, out_wide, sep = "\t", row.names = FALSE,
            quote = FALSE, na = "")
cat(sprintf("Wrote wide table (%d rows): %s\n", nrow(wide), out_wide))

cat("\nSummary:\n")
cat(sprintf("  samples with at least 1 HPO term: %d\n", sum(wide$n_hpo_terms > 0)))
cat(sprintf("  samples with no HPO terms recorded: %d\n", sum(wide$n_hpo_terms == 0)))
cat(sprintf("  median HPO terms per sample: %g\n",
            median(wide$n_hpo_terms[wide$n_hpo_terms > 0])))
