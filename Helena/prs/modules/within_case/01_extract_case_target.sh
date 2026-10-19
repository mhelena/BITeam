#!/usr/bin/env bash
# Cases-only subset of the QC'd target, with MAF/HWE/missingness filters reapplied.
# Usage: bsub < 01_extract_case_target.sh
#BSUB -P re_gecip_neurology
#BSUB -q short -M 8000 -n 2 -R "rusage[mem=8000]"
#BSUB -J within01 -o logs/%J_within01.out -e logs/%J_within01.err
set -euo pipefail
source "${PIPELINE_CONFIG:?export PIPELINE_CONFIG=/path/to/config.sh}"
module load ${PLINK2_MODULE}

IN_PREFIX="${WORKDIR}/merged/case-HC-DEE-EPI.QC"
KEEP="${WORKDIR}/merged/case_samples.txt"
mkdir -p "${WITHIN_DIR}"
OUT_PREFIX="${WITHIN_DIR}/case_only.QC"

echo "input: ${IN_PREFIX}, keep: ${KEEP} ($(wc -l < ${KEEP}) samples)"

plink2 --bfile "${IN_PREFIX}" \
       --keep "${KEEP}" \
       --maf 0.01 --hwe 1e-6 --geno 0.01 --mind 0.01 \
       --make-bed \
       --threads 2 --memory 7000 \
       --out "${OUT_PREFIX}"

echo "cases-only target: $(wc -l < ${OUT_PREFIX}.bim) variants, $(wc -l < ${OUT_PREFIX}.fam) samples"
