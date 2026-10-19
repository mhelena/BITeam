#!/usr/bin/env bash
# PRSice-2 for every score in PRS_LIST, cases vs PRSICE_CONTROL, all samples scored.
# Usage: bsub < 06_run_prsice.sh
#BSUB -P re_gecip_neurology
#BSUB -q long -M 16000 -n 4 -R "rusage[mem=16000]"
#BSUB -J prs06 -o logs/%J_prs06.out -e logs/%J_prs06.err
set -euo pipefail
source "${PIPELINE_CONFIG:?export PIPELINE_CONFIG=/path/to/config.sh}"

cd ${WORKDIR}
mkdir -p results

TARGET=merged/case-HC-DEE-EPI.QC

for phe in ${PRS_LIST}; do
    var=GWAS_${phe}
    GWAS_FILE="${!var:-}"
    [[ -n ${GWAS_FILE} ]] || { echo "no ${var} in config.sh, skipping ${phe}"; continue; }
    [[ -e ${GWAS_FILE} ]] || { echo "MISSING base file: ${GWAS_FILE}"; exit 1; }

    echo "[$(date)] PRSice: ${phe} (base ${GWAS_FILE})"
    "${PRSICE}" \
        --A1 A1 --A2 A2 \
        --all-score \
        --bar-levels 0.0001,0.001,0.01,0.05,0.1,0.5,1 \
        --base "${GWAS_FILE}" \
        --binary-target T \
        --bp BP --chr CHR \
        --fastscore \
        --ignore-fid \
        --model add \
        --pvalue P --se SE \
        --snp SNP --stat BETA \
        --cov-file cov.txt \
        --pheno-file pheno.txt \
        --target ${TARGET} \
        --thread 4 \
        --print-snp \
        --out results/PRS_${phe}
done

echo "PRSice outputs in ${WORKDIR}/results/"
