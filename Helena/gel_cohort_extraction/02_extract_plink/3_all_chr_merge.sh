#!/usr/bin/env bash
#BSUB -q medium
#BSUB -M 12000
#BSUB -n 2
#BSUB -R "rusage[mem=12000]"
# Merges the 22 per-chromosome BEDs, then applies MAF >= 0.05, geno < 0.01 and
# removes A/T and C/G SNPs. Submitted by submit.sh, which sets CHRDIR and FINALDIR.
set -euo pipefail
: "${PIPELINE_CONFIG:?export PIPELINE_CONFIG=/absolute/path/to/config.sh first}"
source "$PIPELINE_CONFIG"

: "${COHORT:?}"
: "${CHRDIR:?}"
: "${FINALDIR:?}"
mkdir -p "$FINALDIR"

module load "$PLINK2_MODULE"

merged="${FINALDIR}/${COHORT}_merged"
mafgeno="${FINALDIR}/${COHORT}_mafgeno"
final="${FINALDIR}/${COHORT}_ALL"
atgc_bed="${FINALDIR}/${COHORT}_atgc.bed"
list="${FINALDIR}/${COHORT}_chr_list.txt"

> "$list"
for c in $(seq 1 22); do
    pref="${CHRDIR}/${COHORT}_chr${c}"
    [[ -s "${pref}.bed" ]] && echo "$pref" >> "$list"
done
n=$(wc -l < "$list" | tr -d ' ')
echo "[$(date)] merging $n per-chr BEDs"

plink2 \
    --pmerge-list "$list" bfile \
    --make-bed \
    --memory 30000 --threads 4 \
    --out "$merged"

plink2 \
    --bfile "$merged" \
    --maf 0.05 \
    --geno 0.01 \
    --make-bed \
    --memory 30000 --threads 4 \
    --out "$mafgeno"

awk 'BEGIN{OFS="\t"} \
     ($5=="A" && $6=="T") || ($5=="T" && $6=="A") || \
     ($5=="C" && $6=="G") || ($5=="G" && $6=="C") {print $1, $4-1, $4}' \
    "${mafgeno}.bim" > "$atgc_bed"

n_atgc=$(wc -l < "$atgc_bed" | tr -d ' ')
echo "[$(date)] A/T and C/G SNPs to remove: $n_atgc"

if [[ "$n_atgc" -gt 0 ]]; then
    plink2 \
        --bfile "$mafgeno" \
        --exclude bed1 "$atgc_bed" \
        --make-bed \
        --memory 30000 --threads 4 \
        --out "$final"
else
    cp "${mafgeno}.bed" "${final}.bed"
    cp "${mafgeno}.bim" "${final}.bim"
    cp "${mafgeno}.fam" "${final}.fam"
fi

rm -f "${merged}".{bed,bim,fam,log} "${mafgeno}".{bed,bim,fam,log}

echo "[$(date)] final: ${final}.{bed,bim,fam}"
wc -l "${final}.bim" "${final}.fam"
