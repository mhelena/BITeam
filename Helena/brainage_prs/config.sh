# Paths read by every script (bash sources it, R and Python parse it).
# Keep one KEY="value" per line with absolute paths; no ${VAR} references.

# Clinical metadata and covariates (columns listed in README.md)
METADATA="/path/to/metadata_minimal.txt"
COVARIATES="/path/to/cov.txt"

# PLINK target genotypes, prefix of .bed/.bim/.fam
TARGET="/path/to/target_QC"

PRSICE="/path/to/PRSice"
# SNP list from an earlier PRSice run (.valid), reused so every base is scored on the same SNPs
VALID_SNPS="/path/to/PRS_epilepsy_syndrome.valid"

# GWAS base files (GRCh37; columns SNP, CHR, BP, A1, A2, BETA, P)
GWAS_epilepsy="/path/to/ILAE3_Epilepsy_37.txt"
GWAS_intelligence="/path/to/Intelligence_37.txt"
GWAS_longevity="/path/to/Longevity_37.txt"
GWAS_AD="/path/to/AD_37.txt"
GWAS_brainage_gm="/path/to/BrainAge_gm.txt"
GWAS_brainage_wm="/path/to/BrainAge_wm.txt"
GWAS_brainage_gwm="/path/to/BrainAge_gwm.txt"
GWAS_frailty="/path/to/Frailty_37_chrpos.txt"

# Outputs: phenotype files (step 0), PRSice scores (step 1), tables and figures (steps 2-7)
PHENO_DIR="/path/to/pheno_files"
PRS_DIR="/path/to/PRS_scores"
OUT_DIR="/path/to/results"
