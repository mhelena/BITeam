#!/usr/bin/env bash
#BSUB -q medium
#BSUB -M 6000
#BSUB -n 2
#BSUB -R "rusage[mem=6000]"
# Merges the chunk BEDs of one chromosome (array index 1-22) with plink2.
# Submitted by submit.sh, which sets OUTDIR and CHRDIR.
set -euo pipefail
: "${PIPELINE_CONFIG:?export PIPELINE_CONFIG=/absolute/path/to/config.sh first}"
source "$PIPELINE_CONFIG"

: "${COHORT:?}"
: "${OUTDIR:?}"
: "${CHRDIR:?}"
mkdir -p "$CHRDIR"

module load "$PLINK2_MODULE"
chr="${LSB_JOBINDEX}"
list="${CHRDIR}/chunks_chr${chr}.txt"

# Chromosome is read from each .bim rather than chunks.txt, so the list does not
# depend on chunk order or naming.
for bim in "${OUTDIR}"/chunk_*.bim; do
    [[ -s "$bim" ]] || continue
    first_chr=$(awk 'NR==1{print $1; exit}' "$bim")
    first_chr="${first_chr#chr}"
    if [[ "$first_chr" == "$chr" ]]; then
        echo "${bim%.bim}"
    fi
done | sort -V > "$list"

n=$(wc -l < "$list" | tr -d ' ')
echo "[$(date)] chr${chr}: merging $n chunks"
if [[ "$n" -eq 0 ]]; then
    echo "[$(date)] chr${chr}: no chunks, nothing to do"
    exit 0
fi
if [[ "$n" -eq 1 ]]; then
    pref=$(cat "$list")
    cp "${pref}.bed" "${CHRDIR}/${COHORT}_chr${chr}.bed"
    cp "${pref}.bim" "${CHRDIR}/${COHORT}_chr${chr}.bim"
    cp "${pref}.fam" "${CHRDIR}/${COHORT}_chr${chr}.fam"
    exit 0
fi

plink2 \
    --pmerge-list "$list" bfile \
    --make-bed \
    --memory 14000 \
    --threads 4 \
    --out "${CHRDIR}/${COHORT}_chr${chr}"

echo "[$(date)] chr${chr}: done"
