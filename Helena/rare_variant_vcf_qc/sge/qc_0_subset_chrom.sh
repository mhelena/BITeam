#!/usr/bin/env bash
#$ -S /bin/bash
#$ -cwd
#$ -j n
#$ -l h_rt=06:00:00
#$ -l h_vmem=4G
#$ -l tmem=4G
# Splits RAW_VCF into one BCF per chromosome (array task = line in CHROM_LIST),
# restricted to SAMPLE_FILE when set. Submitted by qc_submit.sh.
set -euo pipefail

: "${OUTDIR:?}"
: "${RAW_VCF:?}"
: "${CHROM_LIST:?}"
: "${MODE:?}"

PART="${OUTDIR}/_chunks"
mkdir -p "$PART"

chrom=$(sed -n "${SGE_TASK_ID}p" "$CHROM_LIST")
[[ -n "$chrom" ]] || { echo "no chrom at task ${SGE_TASK_ID}"; exit 1; }

out="${PART}/${chrom}.${MODE}.bcf"
done_flag="${PART}/${chrom}.subset.done"

if [[ -e "$done_flag" ]]; then
    echo "[$(date)] ${chrom}: subset already done"; exit 0
fi

if [[ -n "${SAMPLE_FILE:-}" ]]; then
    echo "[$(date)] ${chrom}: subsetting to $(wc -l < "$SAMPLE_FILE" | tr -d ' ') samples"
    "$BCFTOOLS" view --no-version --threads 2 \
            -S "$SAMPLE_FILE" --force-samples \
            -r "$chrom" \
            -Ob -o "${out}.tmp" "$RAW_VCF"
else
    echo "[$(date)] ${chrom}: extracting, all samples"
    "$BCFTOOLS" view --no-version --threads 2 \
            -r "$chrom" \
            -Ob -o "${out}.tmp" "$RAW_VCF"
fi

mv "${out}.tmp" "$out"
"$BCFTOOLS" index -f --threads 2 "$out"

n=$("$BCFTOOLS" index -n "$out")
echo "[$(date)] ${chrom}: ${n} sites"
: > "$done_flag"
