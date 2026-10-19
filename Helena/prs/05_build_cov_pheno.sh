#!/usr/bin/env bash
# Write cov.txt (IID PC1-4 Sex) and pheno.txt (IID Pheno) for PRSice.
# Usage: bash 05_build_cov_pheno.sh
set -euo pipefail
source "${PIPELINE_CONFIG:?export PIPELINE_CONFIG=/path/to/config.sh}"

cd ${WORKDIR}
FAM=merged/case-HC-DEE-EPI.QC.fam
EIG=pca/pca_target_only.eigenvec
CASES=merged/case_samples.txt
ANCHOR=merged/${PRSICE_CONTROL}_samples.txt
[[ -s ${ANCHOR} ]] || { echo "no ${ANCHOR}: check PRSICE_CONTROL in config.sh (epi, hc or dee)"; exit 1; }

# PRSice runs with --ignore-fid, so IID comes first. Sex is taken from
# ANCESTRY_SEX_TSV because the aggV2 genotype files carry 0.
awk -F'\t' -v OFS='\t' '
    FILENAME==ARGV[1] {if (FNR>1) sex[$1]=$4; next}
    FNR==1 {print "IID","PC1","PC2","PC3","PC4","Sex"; next}
    {print $2, $3, $4, $5, $6, (($2 in sex) ? sex[$2] : 0)}' \
    ${ANCESTRY_SEX_TSV} ${EIG} > cov.txt

echo "cov.txt: $(( $(wc -l < cov.txt) - 1 )) samples"
awk 'NR>1 {c[$6]++} END {printf "  sex: male %d, female %d, unknown %d\n", c[1], c[2], c[0]}' cov.txt

# 1 = cases, 0 = PRSICE_CONTROL. The other controls are left out on purpose:
# PRSice still scores them, but they do not influence the threshold it picks.
awk -v OFS='\t' '
    BEGIN {print "IID","Pheno"}
    FILENAME==ARGV[1] {cases[$2]; next}
    FILENAME==ARGV[2] {anchor[$2]; next}
    ($2 in cases)  {print $2, 1; next}
    ($2 in anchor) {print $2, 0}' \
    ${CASES} ${ANCHOR} ${FAM} > pheno.txt

echo "pheno.txt: $(awk 'NR>1 && $2==1' pheno.txt | wc -l | tr -d ' ') cases, $(awk 'NR>1 && $2==0' pheno.txt | wc -l | tr -d ' ') ${PRSICE_CONTROL} controls"
