#!/usr/bin/env bash
# PCs of the QC'd target, used as PRSice covariates.
# Usage: bsub < 04_pca_target_only.sh
#BSUB -P re_gecip_neurology
#BSUB -q short -M 8000 -n 2 -R "rusage[mem=8000]"
#BSUB -J prs04 -o logs/%J_prs04.out -e logs/%J_prs04.err
set -euo pipefail
source "${PIPELINE_CONFIG:?export PIPELINE_CONFIG=/path/to/config.sh}"

module load ${PLINK2_MODULE}
cd ${WORKDIR}
mkdir -p pca

plink2 --bfile merged/case-HC-DEE-EPI.QC \
       --maf 0.01 --geno 0.02 \
       --indep-pairwise 200 50 0.1 \
       --threads 2 --memory 7000 --out pca/prune_target

plink2 --bfile merged/case-HC-DEE-EPI.QC \
       --extract pca/prune_target.prune.in \
       --pca 10 \
       --threads 2 --memory 7000 --out pca/pca_target_only

echo "Wrote pca/pca_target_only.{eigenvec,eigenval}"
