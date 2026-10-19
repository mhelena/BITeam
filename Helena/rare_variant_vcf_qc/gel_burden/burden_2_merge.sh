#!/usr/bin/env bash
#BSUB -q medium
#BSUB -M 12000
#BSUB -n 2
#BSUB -R "rusage[mem=12000]"
# Merges the per-cohort VCFs from burden_1 into ${OUTDIR}/${MERGE_OUT_TAG}_for_burden.vcf.gz.
# Submitted by burden_submit.sh with OUTDIR, COHORTS, MERGE_OUT_TAG.
set -euo pipefail

: "${OUTDIR:?}"
: "${COHORTS:?}"
: "${MERGE_OUT_TAG:?}"

module load bcftools/1.16
module load htslib/1.21

inputs=""
for c in $COHORTS; do
    f="${OUTDIR}/${c}_allchr.vcf.gz"
    [[ -s "$f" ]] || { echo "missing input: $f"; exit 1; }
    inputs="$inputs $f"
done

MERGED="${OUTDIR}/${MERGE_OUT_TAG}_for_burden.vcf.gz"

echo "[$(date)] merging ${COHORTS} -> ${MERGED}"
bcftools merge --threads 2 \
    -m none \
    -Oz -o "$MERGED" \
    $inputs

bcftools index --tbi --threads 2 "$MERGED"

n_var=$(bcftools index --stats "$MERGED" | awk '{s+=$3} END{print s+0}')
n_sam=$(bcftools query -l "$MERGED" | wc -l | tr -d ' ')
echo "[$(date)] merged VCF: ${n_var} variants, ${n_sam} samples"
echo "         path: ${MERGED}"
