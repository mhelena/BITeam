#!/usr/bin/env bash
#$ -S /bin/bash
#$ -cwd
#$ -j n
#$ -l h_rt=06:00:00
#$ -l h_vmem=6G
#$ -l tmem=6G
# Concatenates the per-chromosome QC VCFs, in CHROM_LIST order, into
# ${OUTDIR}/study_qc_${MODE}.vcf.gz.
set -euo pipefail

: "${OUTDIR:?}"
: "${CHROM_LIST:?}"
: "${MODE:?}"

PART="${OUTDIR}/_chunks"
out="${OUTDIR}/study_qc_${MODE}.vcf.gz"
log="${OUTDIR}/study_qc_${MODE}.filters.log"
list="${OUTDIR}/_chunks_to_concat.txt"

: > "$log"
echo "[$(date)] concat start" | tee -a "$log"

: > "$list"
while read -r chrom; do
    f="${PART}/${chrom}.qc.vcf.gz"
    [[ -s "$f" ]] || { echo "missing: $f" | tee -a "$log"; exit 1; }
    echo "$f" >> "$list"
done < "$CHROM_LIST"

n_chunks=$(wc -l < "$list" | tr -d ' ')
echo "[$(date)] concatenating ${n_chunks} per-chrom VCFs" | tee -a "$log"

# --naive is safe here: every chunk comes from the same header and chromosomes do not overlap.
"$BCFTOOLS" concat --threads 4 --naive -f "$list" -Oz -o "$out"
"$TABIX" -p vcf -f "$out"

n_vars=$("$BCFTOOLS" index -n "$out")
n_samples=$("$BCFTOOLS" query -l "$out" | wc -l | tr -d ' ')
echo "[$(date)] wrote ${out}" | tee -a "$log"
echo "         variants: ${n_vars}"   | tee -a "$log"
echo "         samples : ${n_samples}" | tee -a "$log"

echo "[$(date)] chromosome coverage:" | tee -a "$log"
"$TABIX" -l "$out" | tee -a "$log"
