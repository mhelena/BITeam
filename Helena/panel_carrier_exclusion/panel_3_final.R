#!/usr/bin/env Rscript
# Optional: when POOL_FILE is the union of several lists, report how many of each
# list survive the carrier filter and write the final list as a PLINK --keep file.
# Usage: source "$PIPELINE_CONFIG" && Rscript panel_3_final.R
suppressPackageStartupMessages(library(Rlabkey))

env <- function(k) {
  v <- Sys.getenv(k)
  if (!nzchar(v)) stop(k, " is not set; source config.sh first")
  v
}
SUBSETS  <- strsplit(trimws(env("POOL_SUBSETS")), "[[:space:]]+")[[1]]
CARRIERS <- file.path(env("WORK_DIR"), "panel_carriers.txt")
OUT_DIR  <- env("FINAL_DIR")
TAG      <- sprintf("panel%s_v%s", env("PANEL_ID"), env("PANEL_VERSION"))
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

stopifnot(all(file.exists(SUBSETS)), file.exists(CARRIERS))

# Lists may hold participant ids rather than platekeys; map through LabKey.
anc <- labkey.selectRows(baseUrl = env("LABKEY_BASE_URL"), folderPath = env("LABKEY_FOLDER"),
                         schemaName = "lists", queryName = "aggregate_gvcf_sample_stats",
                         colNameOpt = "rname", maxRows = 1e7)
PK        <- intersect(c("platekey","plate_key","sample_id"), names(anc))[1]
pid_set   <- unique(as.character(anc$participant_id))
pid_to_pk <- setNames(as.character(anc[[PK]]), as.character(anc$participant_id))

as_platekeys <- function(f) {
  x <- trimws(readLines(f)); x <- x[nzchar(x)]
  if (mean(x %in% pid_set) > 0.5) unname(pid_to_pk[x]) else x
}

pk_carriers <- unique(trimws(readLines(CARRIERS))); pk_carriers <- pk_carriers[nzchar(pk_carriers)]

write_pk <- function(v, f) {
  writeLines(as.character(v), f)
  cat(sprintf("  wrote %d platekeys -> %s\n", length(v), f))
}

cat(sprintf("\nSurvival after %s carrier filter (%d carriers)\n", TAG, length(pk_carriers)))
pk_final <- character(0)
for (f in SUBSETS) {
  pk_in   <- unique(na.omit(as_platekeys(f)))
  pk_kept <- setdiff(pk_in, pk_carriers)
  # union in list order, so earlier lists keep their position
  pk_final <- union(pk_final, pk_kept)
  name <- sub("\\.[^.]*$", "", basename(f))
  cat(sprintf("  %-40s kept %d / %d (%.1f%%)\n", name,
              length(pk_kept), length(pk_in),
              100 * length(pk_kept) / max(1, length(pk_in))))
  write_pk(pk_kept, file.path(OUT_DIR, sprintf("%s_kept.platekeys.txt", name)))
}
cat(sprintf("  final controls (union, carrier-filtered): %d\n", length(pk_final)))

write_pk(pk_final, file.path(OUT_DIR, sprintf("controls_%s.platekeys.txt", TAG)))
write.table(data.frame(FID = pk_final, IID = pk_final),
            file.path(OUT_DIR, sprintf("controls_%s.plink.keep", TAG)),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
cat(sprintf("  wrote PLINK --keep -> %s/controls_%s.plink.keep\n", OUT_DIR, TAG))
