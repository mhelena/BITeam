#!/usr/bin/env bash
#BSUB -q medium
#BSUB -M 16000
#BSUB -n 2
#BSUB -R "rusage[mem=16000]"
# Merges one cohort's 22 per-chromosome BEDs (no MAF filter, rare variants kept)
# and exports a bgzipped VCF. Submitted by burden_submit.sh with COHORT, PERCHR, OUTDIR.
set -euo pipefail

: "${COHORT:?}"
: "${PERCHR:?}"
: "${OUTDIR:?}"
mkdir -p "$OUTDIR"

module load plink/2.0
module load bcftools/1.16
module load htslib/1.21

WORK="${OUTDIR}/_work_${COHORT}"
mkdir -p "$WORK"

LIST="${WORK}/${COHORT}_chrs.txt"
: > "$LIST"
for chr in $(seq 1 22); do
    pref="${PERCHR}/${COHORT}_chr${chr}"
    if [[ -s "${pref}.bed" ]]; then
        echo "$pref" >> "$LIST"
    fi
done
n=$(wc -l < "$LIST" | tr -d ' ')
echo "[$(date)] ${COHORT}: merging $n per-chr BEDs"

plink2 --pmerge-list "$LIST" bfile \
       --make-bed \
       --threads 2 --memory 14000 \
       --out "${WORK}/${COHORT}_allchr"

n_var=$(wc -l < "${WORK}/${COHORT}_allchr.bim" | tr -d ' ')
n_sam=$(wc -l < "${WORK}/${COHORT}_allchr.fam" | tr -d ' ')
echo "[$(date)] ${COHORT}: merged BED has ${n_var} variants, ${n_sam} samples"

plink2 --bfile "${WORK}/${COHORT}_allchr" \
       --export vcf-4.2 bgz id-paste=iid \
       --threads 2 --memory 14000 \
       --out "${OUTDIR}/${COHORT}_allchr"

bcftools index --tbi --threads 2 "${OUTDIR}/${COHORT}_allchr.vcf.gz"

echo "[$(date)] ${COHORT}: wrote ${OUTDIR}/${COHORT}_allchr.vcf.gz"
# The merged cohort BED in ${WORK} is kept; delete it if disk is tight.
