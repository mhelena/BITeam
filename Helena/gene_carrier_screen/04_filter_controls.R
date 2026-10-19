#!/usr/bin/env Rscript
# Keeps carriers from the whitelisted disease categories, optionally EUR >= 0.8, and writes four lists.
# Usage: source $PIPELINE_CONFIG; Rscript 04_filter_controls.R [id_list]   (default: $CARRIERS)
suppressPackageStartupMessages(library(Rlabkey))

cfg <- function(x) {
  v <- Sys.getenv(x)
  if (!nzchar(v)) stop(x, " is not set; source config.sh first", call. = FALSE)
  v
}

args      <- commandArgs(trailingOnly = TRUE)
IN_FILE   <- if (length(args) >= 1) args[1] else cfg("CARRIERS")
WHITELIST <- cfg("WHITELIST")
EXCLUDE   <- Sys.getenv("EXCLUDE_IDS")
OUT_DIR   <- dirname(IN_FILE)
TAG       <- tools::file_path_sans_ext(basename(IN_FILE))

stopifnot("input file not found"     = file.exists(IN_FILE))
stopifnot("whitelist file not found" = file.exists(WHITELIST))
if (nzchar(EXCLUDE)) stopifnot("exclusion file not found" = file.exists(EXCLUDE))

# Input may have a header or several columns; the ID column is chosen after LabKey loads.
raw_lines <- trimws(readLines(IN_FILE))
raw_lines <- raw_lines[nchar(raw_lines) > 0]
cat(sprintf("[%s] input rows: %d\n", format(Sys.time()), length(raw_lines)))
cat(sprintf("  first line: '%s'\n", raw_lines[1]))

parts  <- strsplit(raw_lines, "\\s+")
n_cols <- max(sapply(parts, length))
cat(sprintf("  max columns: %d\n", n_cols))

fetch <- function(q) {
  cat(sprintf("  loading %s\n", q))
  labkey.selectRows(baseUrl = cfg("LABKEY_BASE"), folderPath = cfg("LABKEY_FOLDER"),
                    schemaName = "lists", queryName = q,
                    colNameOpt = "rname", maxRows = 1e7)
}
anc  <- fetch("aggregate_gvcf_sample_stats")
dz   <- fetch("rare_diseases_participant_disease")
hpo  <- fetch("rare_diseases_participant_phenotype")
hpo_present <- hpo[hpo$hpo_present %in% c(TRUE,"Yes","yes","Y",1,"1"), ]

PK <- intersect(c("platekey","plate_key","sample_id"), names(anc))[1]
stopifnot(!is.na(PK))
pk_set  <- unique(as.character(anc[[PK]]))
pid_set <- unique(as.character(anc$participant_id))

best <- list(col = NA, kind = NA, hits = 0, ids = character(0))
for (col in seq_len(n_cols)) {
  ids_col <- unique(na.omit(sapply(parts, function(x)
                                   if (length(x) >= col) x[col] else NA)))
  h_pk  <- sum(ids_col %in% pk_set)
  h_pid <- sum(ids_col %in% pid_set)
  cat(sprintf("  col %d  ->  matches: platekey=%d  participant_id=%d  (n_ids=%d)\n",
              col, h_pk, h_pid, length(ids_col)))
  if (h_pk >= max(best$hits, h_pid)) {
    best <- list(col = col, kind = "platekey",       hits = h_pk,  ids = ids_col)
  }
  if (h_pid > best$hits) {
    best <- list(col = col, kind = "participant_id", hits = h_pid, ids = ids_col)
  }
}

if (best$hits == 0) {
  cat("\nNo matches in any column. Check the input file format.\n")
  cat("First 5 input lines:\n"); print(head(raw_lines, 5))
  cat("First 5 LabKey platekeys:\n"); print(head(pk_set, 5))
  cat("First 5 LabKey participant_ids:\n"); print(head(pid_set, 5))
  quit(status = 1)
}

cat(sprintf("\n  -> using column %d as %s (%d matches)\n",
            best$col, best$kind, best$hits))
ids_in <- best$ids

pk_to_pid <- setNames(as.character(anc$participant_id), as.character(anc[[PK]]))
if (best$kind == "participant_id") {
  pids_in <- intersect(ids_in, pid_set)
} else {
  pids_in <- na.omit(unique(pk_to_pid[ids_in]))
}
cat(sprintf("  participant_ids in input: %d\n", length(pids_in)))

if (nzchar(EXCLUDE)) {
  excl <- trimws(readLines(EXCLUDE))
  excl <- excl[nchar(excl) > 0]
  excl_pids <- union(intersect(excl, pid_set), na.omit(unname(pk_to_pid[excl])))
  n0 <- length(pids_in)
  pids_in <- setdiff(pids_in, excl_pids)
  cat(sprintf("  exclusion list: %d IDs, %d removed from input\n",
              length(excl), n0 - length(pids_in)))
}

terms <- trimws(readLines(WHITELIST))
terms <- terms[nchar(terms) > 0]
lo <- tolower(terms)
dz_hits <- unique(dz$participant_id[
  tolower(dz$specific_disease)  %in% lo |
  tolower(dz$disease_sub_group) %in% lo |
  tolower(dz$disease_group)     %in% lo])
hpo_hits <- unique(hpo_present$participant_id[
  tolower(hpo_present$normalised_hpo_term) %in% lo])
whitelist_ids <- union(dz_hits, hpo_hits)
cat(sprintf("\n[whitelist] terms: %d   matched participants in GEL: %d\n",
            length(terms), length(whitelist_ids)))

# The list holds categories unrelated to the phenotype, so matches are kept (intersect,
# not setdiff).
pids_wl   <- intersect(pids_in, whitelist_ids)
pids_drop <- setdiff(pids_in, whitelist_ids)
cat(sprintf("[whitelist] kept: %d   dropped (no match): %d\n",
            length(pids_wl), length(pids_drop)))

pid_to_pk <- setNames(as.character(anc[[PK]]), as.character(anc$participant_id))
pids_to_pks <- function(pids) {
  pks <- pid_to_pk[as.character(pids)]
  unique(na.omit(pks))
}

prob_eur <- intersect(c("pred_european_ancestries","pred_european_prob",
                        "european_prob"), names(anc))[1]
stopifnot("No EUR ancestry column in aggregate_gvcf_sample_stats" = !is.na(prob_eur))
pid_to_eur <- setNames(as.numeric(anc[[prob_eur]]), as.character(anc$participant_id))

eur_filter <- function(pids) {
  v <- pid_to_eur[as.character(pids)]
  pids[!is.na(v) & v >= 0.8]
}
pids_eur08    <- eur_filter(pids_in)
pids_wl_eur08 <- eur_filter(pids_wl)
cat(sprintf("\n[EUR>=0.8] input          : %d  -> %d kept\n",
            length(pids_in), length(pids_eur08)))
cat(sprintf("[EUR>=0.8] after whitelist: %d  -> %d kept\n",
            length(pids_wl), length(pids_wl_eur08)))

out_allanc   <- file.path(OUT_DIR, sprintf("%s_allanc.txt",          TAG))
out_wl       <- file.path(OUT_DIR, sprintf("%s_whitelist.txt",       TAG))
out_eur08    <- file.path(OUT_DIR, sprintf("%s_eur08.txt",           TAG))
out_wl_eur08 <- file.path(OUT_DIR, sprintf("%s_whitelist_eur08.txt", TAG))

pks_in <- pids_to_pks(pids_in)
writeLines(pks_in,                     out_allanc)
writeLines(pids_to_pks(pids_wl),       out_wl)
writeLines(pids_to_pks(pids_eur08),    out_eur08)
writeLines(pids_to_pks(pids_wl_eur08), out_wl_eur08)

cat(sprintf("\nOutputs (in %s):\n", OUT_DIR))
cat(sprintf("  %s   (n=%d, no filters)\n",              basename(out_allanc),   length(pks_in)))
cat(sprintf("  %s   (n=%d, whitelist only)\n",          basename(out_wl),       length(pids_wl)))
cat(sprintf("  %s   (n=%d, EUR>=0.8 only)\n",           basename(out_eur08),    length(pids_eur08)))
cat(sprintf("  %s   (n=%d, whitelist + EUR>=0.8)\n",    basename(out_wl_eur08), length(pids_wl_eur08)))
