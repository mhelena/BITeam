#!/usr/bin/env bash
# PCs within the cases only, same recipe as 04_pca_target_only.sh.
# Usage: bsub < 02_case_pca.sh
#BSUB -P re_gecip_neurology
#BSUB -q short -M 8000 -n 2 -R "rusage[mem=8000]"
#BSUB -J within02 -o logs/%J_within02.out -e logs/%J_within02.err
set -euo pipefail
source "${PIPELINE_CONFIG:?export PIPELINE_CONFIG=/path/to/config.sh}"
module load ${PLINK2_MODULE}

TARGET="${WITHIN_DIR}/case_only.QC"

plink2 --bfile "${TARGET}" \
       --maf 0.01 --geno 0.02 \
       --indep-pairwise 200 50 0.1 \
       --threads 2 --memory 7000 \
       --out "${WITHIN_DIR}/_prune_cases"

plink2 --bfile "${TARGET}" \
       --extract "${WITHIN_DIR}/_prune_cases.prune.in" \
       --pca 10 \
       --threads 2 --memory 7000 \
       --out "${WITHIN_DIR}/pca_case_only"

echo "Wrote ${WITHIN_DIR}/pca_case_only.{eigenvec,eigenval}"
