#!/usr/bin/env Rscript
# Candidate pool: hg38 probands without the phenotype (any ancestry), as platekeys.
# Usage: source $PIPELINE_CONFIG; Rscript 01_build_pool.R   (login node, Rlabkey via ~/.netrc)
suppressPackageStartupMessages(library(Rlabkey))

cfg <- function(x) {
  v <- Sys.getenv(x)
  if (!nzchar(v)) stop(x, " is not set; source config.sh first", call. = FALSE)
  v
}

SAMPLES <- cfg("POOL_FILE")
dir.create(dirname(SAMPLES), showWarnings = FALSE, recursive = TRUE)

fetch <- function(q) {
  cat(sprintf("[%s] loading %s\n", format(Sys.time()), q))
  labkey.selectRows(baseUrl = cfg("LABKEY_BASE"), folderPath = cfg("LABKEY_FOLDER"),
                    schemaName = "lists", queryName = q,
                    colNameOpt = "rname", maxRows = 1e7)
}

anc    <- fetch("aggregate_gvcf_sample_stats")
psum   <- fetch("participant_summary")
dz     <- fetch("rare_diseases_participant_disease")
part   <- fetch("participant")
hpo    <- fetch("rare_diseases_participant_phenotype")
gfp    <- fetch("genome_file_paths_and_types")

hg38_ids    <- unique(gfp$participant_id[
                 grepl("38|GRCh38", gfp$genome_build, ignore.case = TRUE)])
proband_ids <- unique(part$participant_id[part$participant_type == "Proband"])

pheno_regex <- cfg("PHENO_HPO_REGEX")
pheno_sub   <- trimws(strsplit(Sys.getenv("PHENO_DISEASE_SUBGROUP"), ";")[[1]])
pheno_sub   <- unique(c(pheno_sub, tolower(pheno_sub)))
pheno_sub   <- pheno_sub[nzchar(pheno_sub)]
hpo_present <- hpo[hpo$hpo_present %in% c(TRUE,"Yes","yes","Y",1,"1"), ]
pheno_disease <- unique(dz$participant_id[dz$disease_sub_group %in% pheno_sub])
hpo_pheno <- unique(hpo_present$participant_id[
              grepl(pheno_regex, hpo_present$normalised_hpo_term, ignore.case = TRUE)])
psum_pheno <- unique(psum$participant_id[
               grepl(pheno_regex, psum$hpo_term, ignore.case = TRUE)])
pheno_ids <- Reduce(union, list(pheno_disease, hpo_pheno, psum_pheno))

cat(sprintf("\nProband participants (hg38): %d\n",
            length(intersect(proband_ids, hg38_ids))))
cat(sprintf("Participants with the phenotype: %d\n", length(pheno_ids)))

candidate_ids <- setdiff(intersect(proband_ids, hg38_ids), pheno_ids)
cat(sprintf("Probands on hg38 without the phenotype: %d\n", length(candidate_ids)))

PK <- intersect(c("platekey","plate_key","sample_id"), names(anc))[1]
stopifnot("No platekey column in aggregate_gvcf_sample_stats" = !is.na(PK))
pk_map <- unique(anc[, c("participant_id", PK)])
pk_map <- pk_map[!is.na(pk_map$participant_id) & !is.na(pk_map[[PK]]), ]
pks <- pk_map[[PK]][match(as.character(candidate_ids),
                          as.character(pk_map$participant_id))]
pks <- unique(na.omit(pks))
cat(sprintf("Mapped to platekeys (pool size): %d\n", length(pks)))

# Ancestry is reported here but only filtered later, in 04_filter_controls.R.
prob_eur <- intersect(c("pred_european_ancestries","pred_european_prob",
                        "european_prob"), names(anc))[1]
if (!is.na(prob_eur)) {
  pk_to_eur <- setNames(as.numeric(anc[[prob_eur]]), as.character(anc[[PK]]))
  eur_vals  <- pk_to_eur[pks]
  report <- function(label, keep) {
    cat(sprintf("  %-20s : %d  (%.1f%%)\n", label, sum(keep), 100 * sum(keep) / length(pks)))
  }
  cat("\nAncestry breakdown (not filtered):\n")
  report("EUR >= 0.8",       !is.na(eur_vals) & eur_vals >= 0.8)
  report("EUR 0.5 - 0.8",    !is.na(eur_vals) & eur_vals >= 0.5 & eur_vals < 0.8)
  report("EUR < 0.5",        !is.na(eur_vals) & eur_vals < 0.5)
  report("EUR prob missing", is.na(eur_vals))
}

writeLines(pks, SAMPLES)
cat(sprintf("\nWrote %s (%d platekeys).\n", SAMPLES, length(pks)))
