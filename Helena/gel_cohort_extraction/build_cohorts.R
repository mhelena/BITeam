#!/usr/bin/env Rscript
# Builds the EUR >= 0.8 control cohorts from LabKey and writes <out_dir>/<cohort>/samples.txt (platekeys).
# Usage: Rscript build_cohorts.R <out_dir> [non_epi_terms.txt]    (needs Rlabkey and a LabKey key in ~/.netrc)
suppressPackageStartupMessages(library(Rlabkey))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("usage: Rscript build_cohorts.R <out_dir> [non_epi_terms.txt]")
out_root     <- args[1]
non_epi_path <- if (length(args) >= 2) args[2] else NA

base   <- Sys.getenv("LABKEY_BASE", "https://labkey.prod.aws.gel.ac/labkey")
folder <- Sys.getenv("LABKEY_FOLDER", "/main-programme/main-programme_v19_2024-10-31")
schema <- "lists"

cohorts <- c("epilepsy_controls", "dee", "healthy_pool_candidate", "scn1a_pool")
for (d in cohorts) dir.create(file.path(out_root, d), showWarnings = FALSE, recursive = TRUE)

fetch <- function(name) {
  cat(sprintf("[%s] loading %s ...\n", format(Sys.time()), name))
  labkey.selectRows(baseUrl = base, folderPath = folder,
                    schemaName = schema, queryName = name,
                    colNameOpt = "rname", maxRows = 1e7)
}

anc    <- fetch("aggregate_gvcf_sample_stats")
psum   <- fetch("participant_summary")
dz     <- fetch("rare_diseases_participant_disease")
part   <- fetch("participant")
hpo    <- fetch("rare_diseases_participant_phenotype")
solved <- fetch("gmc_exit_questionnaire")
gfp    <- fetch("genome_file_paths_and_types")

eur08    <- unique(anc$participant_id[anc$pred_european_ancestries >= 0.8])
hg38_ids <- unique(gfp$participant_id[grepl("38|GRCh38", gfp$genome_build, ignore.case = TRUE)])

proband_ids    <- unique(part$participant_id[part$participant_type == "Proband"])
relative_ids   <- unique(part$participant_id[part$participant_type == "Relative"])
unaffected_ids <- unique(psum$participant_id[
  psum$affection_status %in% c("Unaffected", "unaffected", "UNAFFECTED")])

solved_ids <- unique(solved$participant_id[
  tolower(as.character(solved$case_solved_family)) == "yes"])

EPI_SUB  <- c("Inherited Epilepsy Syndromes", "inherited epilepsy syndromes")
DEE_SPEC <- c("Epileptic encephalopathy", "epileptic encephalopathy")

epi_disease <- unique(dz$participant_id[dz$disease_sub_group %in% EPI_SUB])
ee_ids      <- unique(dz$participant_id[dz$specific_disease %in% DEE_SPEC])

hpo_present <- hpo[hpo$hpo_present %in% c(TRUE, "Yes", "yes", "Y", 1, "1"), ]
hpo_epi <- unique(hpo_present$participant_id[
  grepl("epilep|seizure|convuls", hpo_present$normalised_hpo_term, ignore.case = TRUE)])
psum_epi <- unique(psum$participant_id[
  grepl("epilep|seizure|convuls", psum$hpo_term, ignore.case = TRUE)])

# Anyone with an epilepsy disease label or an epilepsy/seizure HPO term.
mega_pool <- Reduce(union, list(epi_disease, hpo_epi, psum_epi))

epi_ctrl <- Reduce(intersect, list(
  setdiff(setdiff(mega_pool, ee_ids), solved_ids),
  proband_ids, eur08, hg38_ids))

dee <- Reduce(intersect, list(ee_ids, proband_ids, eur08, hg38_ids))

scn1a_pool <- setdiff(
  Reduce(intersect, list(proband_ids, eur08, hg38_ids)),
  mega_pool)

fam_lookup <- setNames(as.character(part$rare_diseases_family_id),
                       as.character(part$participant_id))

restrict_one_per_fam <- function(ids) {
  df <- data.frame(pid = as.character(ids),
                   fam = fam_lookup[as.character(ids)],
                   stringsAsFactors = FALSE)
  df <- df[!is.na(df$fam) & df$fam != "", ]
  df <- df[order(df$fam, df$pid), ]
  df$pid[!duplicated(df$fam)]
}

unrel_ids <- unique(anc$participant_id[
  anc$set_of_unrelated_individuals %in% c(1, "1", TRUE, "TRUE", "Y", "yes")])

# Healthy candidates: unaffected rare-disease relatives outside the epilepsy pool,
# in GEL's unrelated set, one per family. The final healthy control set also drops
# carriers of gnomAD-absent variants in the PanelApp "Early onset or syndromic
# epilepsy" panel; that filter is a separate pipeline.
healthy_base <- Reduce(intersect, list(
  setdiff(intersect(relative_ids, unaffected_ids), mega_pool), eur08, hg38_ids))
healthy_candidate <- restrict_one_per_fam(intersect(healthy_base, unrel_ids))
cat(sprintf("[healthy candidates] %d\n", length(healthy_candidate)))

# Optional: drop candidates with any disease or HPO term from Supplementary Table 3
# of PMID:37006128 (conditions unrelated to epilepsy), so the pool is free of
# documented medical conditions. One term per line.
if (!is.na(non_epi_path)) {
  if (!file.exists(non_epi_path)) stop("not found: ", non_epi_path)
  terms <- trimws(readLines(non_epi_path))
  terms <- terms[nchar(terms) > 0]
  lo <- tolower(terms)

  dz_hits <- unique(dz$participant_id[
    tolower(dz$specific_disease)  %in% lo |
    tolower(dz$disease_sub_group) %in% lo |
    tolower(dz$disease_group)     %in% lo
  ])
  hpo_hits <- unique(hpo_present$participant_id[
    tolower(hpo_present$normalised_hpo_term) %in% lo
  ])
  supp3_ids <- union(dz_hits, hpo_hits)

  before <- length(healthy_candidate)
  healthy_candidate <- setdiff(healthy_candidate, supp3_ids)
  cat(sprintf("[non-epilepsy term exclusion] terms=%d matched=%d removed=%d remaining=%d\n",
              length(terms), length(supp3_ids), before - length(healthy_candidate),
              length(healthy_candidate)))
} else {
  cat("[non-epilepsy term exclusion] no term file given, skipped\n")
}

cat("\nCohort sizes (participants, EUR >= 0.8, GRCh38)\n")
cat(sprintf("  epilepsy_controls      : %d\n", length(epi_ctrl)))
cat(sprintf("  dee                    : %d\n", length(dee)))
cat(sprintf("  healthy_pool_candidate : %d\n", length(healthy_candidate)))
cat(sprintf("  scn1a_pool             : %d\n", length(scn1a_pool)))

PK_COL <- intersect(c("platekey", "plate_key", "sample_id"), names(anc))[1]
if (is.na(PK_COL))
  stop("No platekey column in aggregate_gvcf_sample_stats, check names(anc)")

map <- unique(anc[, c("participant_id", PK_COL)])
map <- map[!is.na(map[[PK_COL]]) & !is.na(map$participant_id), ]

to_platekeys <- function(pids) {
  pks <- map[[PK_COL]][match(as.character(pids), as.character(map$participant_id))]
  na.omit(pks)
}

write_cohort <- function(nm, pids) {
  pks  <- to_platekeys(pids)
  miss <- length(pids) - length(pks)
  out  <- file.path(out_root, nm, "samples.txt")
  writeLines(as.character(pks), out)
  cat(sprintf("  %-22s participants=%d platekeys=%d missing=%d -> %s\n",
              nm, length(pids), length(pks), miss, out))
}

cat("\nWriting samples.txt\n")
write_cohort("epilepsy_controls",      epi_ctrl)
write_cohort("dee",                    dee)
write_cohort("healthy_pool_candidate", healthy_candidate)
write_cohort("scn1a_pool",             scn1a_pool)

saveRDS(list(epi_ctrl = epi_ctrl, dee = dee, healthy_candidate = healthy_candidate,
             scn1a_pool = scn1a_pool, eur08 = eur08, hg38_ids = hg38_ids,
             proband_ids = proband_ids, mega_pool = mega_pool, ee_ids = ee_ids,
             solved_ids = solved_ids, map = map),
        file.path(out_root, "cohort_objects.rds"))
cat("\nSaved", file.path(out_root, "cohort_objects.rds"), "\n")
