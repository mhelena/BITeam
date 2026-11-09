# Settings shared by every notebook.
# Only the first block normally needs editing (where plink2 and PRSice live).

# 1. Tools (edit these)
# Full path to each program in quotes, or just its name if it is on your PATH,
# e.g. "/Users/me/tools/plink2" or "C:/tools/plink2.exe" (forward slashes on Windows).
PLINK2 <- "plink2"
PRSICE <- "PRSice_linux"     # PRSice_mac on a Mac, PRSice_win64.exe on Windows
THREADS <- 2

# 2. Input data
# For the workshop these point at the toy cohort made by 00_make_toy_data.R.
# For a real cohort, point them at your own files (same formats).
DATA_DIR <- "data/toy"
GENO     <- file.path(DATA_DIR, "toy_cohort")             # PLINK prefix: .bed/.bim/.fam
CLINICAL <- file.path(DATA_DIR, "clinical.tsv")           # one row per person
FINDINGS <- file.path(DATA_DIR, "genetic_findings.tsv")   # one row per reported variant
GWAS <- c(                                                # base GWAS summary statistics
  epilepsy = file.path(DATA_DIR, "gwas_epilepsy.tsv.gz"),
  height   = file.path(DATA_DIR, "gwas_height.tsv.gz")    # negative control: no link to epilepsy
)

# 3. Gene groups (used to define genotype groups)
GENES_HIGH_PENETRANCE <- c("SCN1A", "SCN2A", "KCNQ2", "STXBP1", "CDKL5")
GENES_INCOMPLETE      <- c("DEPDC5", "NPRL2", "NPRL3")

# 4. PRS settings
PT_LEVELS <- c(5e-8, 1e-5, 1e-3, 0.01, 0.05, 0.1, 0.5, 1)   # P value thresholds scored
PT_MAIN   <- 0.1      # the ONE threshold for the main analysis, chosen before looking
N_PCS     <- 4        # principal components used as covariates

# 5. Output folders
RESULTS <- "results"
QC_DIR  <- file.path(RESULTS, "qc")
PRS_DIR <- file.path(RESULTS, "prs")
for (d in c(RESULTS, QC_DIR, PRS_DIR)) dir.create(d, showWarnings = FALSE, recursive = TRUE)
