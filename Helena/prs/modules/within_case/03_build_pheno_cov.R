#!/usr/bin/env Rscript
# PRSice phenotype and covariate files for the within-case analyses (binary outcome, optional onset).
# Usage: Rscript 03_build_pheno_cov.R   (settings from config.sh via the environment, or the defaults below)
suppressPackageStartupMessages({
    library(dplyr)
})

WORKDIR    <- Sys.getenv("WORKDIR", "/re_gecip/<your_gecip>/<project>/<cohort>")
WORK       <- Sys.getenv("WITHIN_DIR", file.path(WORKDIR, "within_case"))
OUTCOME    <- Sys.getenv("OUTCOME", "outcome")
OUTCOME_IN <- Sys.getenv("OUTCOME_TSV", file.path(WORKDIR, "phenotypes/outcome.tsv"))
ONSET_IN   <- Sys.getenv("ONSET_TSV", "")
SEX_IN     <- Sys.getenv("WITHIN_SEX_TSV", file.path(WORKDIR, "phenotypes/sex.tsv"))
EIG_IN     <- file.path(WORK, "pca_case_only.eigenvec")
FAM_IN     <- file.path(WORK, "case_only.QC.fam")

for (f in c(OUTCOME_IN, SEX_IN, EIG_IN, FAM_IN, if (nzchar(ONSET_IN)) ONSET_IN))
    if (!file.exists(f)) stop(sprintf("MISSING: %s", f))

sex_df <- read.table(SEX_IN, header = TRUE, sep = "", stringsAsFactors = FALSE)
sex_df$Sex <- ifelse(tolower(as.character(sex_df$Sex)) %in% c("male","m","1"), 1L,
              ifelse(tolower(as.character(sex_df$Sex)) %in% c("female","f","2"), 2L, NA_integer_))
sex_df <- sex_df[, c("IID","Sex")]

eig <- read.table(EIG_IN, header = TRUE, comment.char = "", check.names = FALSE)
names(eig) <- gsub("^#", "", names(eig))
names(eig)[tolower(names(eig)) == "iid"] <- "IID"
eig_use <- eig[, c("IID","PC1","PC2","PC3","PC4")]

cov <- eig_use %>%
    left_join(sex_df, by = "IID") %>%
    mutate(Sex = ifelse(is.na(Sex), 0L, Sex))
write.table(cov, file.path(WORK, "cov.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("cov.txt: %d rows (M=%d, F=%d, unk=%d)\n",
            nrow(cov),
            sum(cov$Sex == 1), sum(cov$Sex == 2), sum(cov$Sex == 0)))

# PRSice binary coding: yes = 2 (case), no = 1 (control), anything else dropped
oc <- read.table(OUTCOME_IN, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
col <- names(oc)[tolower(names(oc)) == tolower(OUTCOME)][1]
if (is.na(col)) stop(sprintf("no column named %s in %s", OUTCOME, OUTCOME_IN))
oc$Outcome <- ifelse(tolower(oc[[col]]) == "yes", 2L,
              ifelse(tolower(oc[[col]]) == "no",  1L, NA_integer_))
oc_out <- oc %>% filter(!is.na(Outcome)) %>% select(IID, Outcome)
write.table(oc_out, file.path(WORK, "pheno_binary.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("pheno_binary.txt (%s): %d rows (yes=%d, no=%d)\n",
            OUTCOME, nrow(oc_out), sum(oc_out$Outcome == 2), sum(oc_out$Outcome == 1)))

fam <- read.table(FAM_IN, header = FALSE, stringsAsFactors = FALSE)
fam_iids <- fam[[2]]
cat(sprintf("\nSample intersection with the cases-only target (n=%d):\n", length(fam_iids)))
cat(sprintf("  %s covered: %d / %d\n", OUTCOME, sum(oc_out$IID %in% fam_iids), nrow(oc_out)))

if (nzchar(ONSET_IN)) {
    onset <- read.table(ONSET_IN, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    onset$onset_months <- suppressWarnings(as.numeric(onset$onset_months))
    onset_out <- onset %>% filter(!is.na(onset_months)) %>% select(IID, onset_months)
    write.table(onset_out, file.path(WORK, "pheno_onset.txt"),
                sep = "\t", row.names = FALSE, quote = FALSE)
    cat(sprintf("pheno_onset.txt: %d rows (range %.1f - %.1f months, median %.1f)\n",
                nrow(onset_out),
                min(onset_out$onset_months), max(onset_out$onset_months),
                median(onset_out$onset_months)))
    cat(sprintf("  onset covered: %d / %d\n",
                sum(onset_out$IID %in% fam_iids), nrow(onset_out)))
}
