#!/usr/bin/env bash
#$ -S /bin/bash
#$ -cwd
#$ -j n
#$ -l h_rt=06:00:00
#$ -l h_vmem=4G
#$ -l tmem=4G
# Writes the filtered VCF for one chromosome: drops union fail sites and indels
# longer than 50 bp, and LCR_BED intervals when set. No MAF, LD or allele-count filter.
set -euo pipefail

: "${OUTDIR:?}"
: "${CHROM_LIST:?}"
: "${MODE:?}"

PART="${OUTDIR}/_chunks"
chrom=$(sed -n "${SGE_TASK_ID}p" "$CHROM_LIST")
[[ -n "$chrom" ]] || { echo "no chrom at task ${SGE_TASK_ID}"; exit 1; }

in_bcf="${PART}/${chrom}.${MODE}.bcf"
regions="${OUTDIR}/all_fail_regions.tsv.gz"
out="${PART}/${chrom}.qc.vcf.gz"
done_flag="${PART}/${chrom}.filter.done"

[[ -s "$in_bcf"  ]] || { echo "missing input: $in_bcf";  exit 1; }
[[ -s "$regions" ]] || { echo "missing regions: $regions"; exit 1; }

if [[ -e "$done_flag" ]]; then
    echo "[$(date)] ${chrom}: filter already done"; exit 0
fi

if [[ -n "${LCR_BED:-}" ]]; then
    [[ -s "$LCR_BED" ]] || { echo "missing LCR bed: $LCR_BED"; exit 1; }
    # bcftools takes only one -T, so the LCR mask goes in a second view.
    echo "[$(date)] ${chrom}: filtering (exclude fails + LCR + indel<=50bp)"
    "$BCFTOOLS" view --no-version --threads 2 \
            -T "^${regions}" \
            -e 'abs(strlen(REF)-strlen(ALT))>50' \
            -Ou "$in_bcf" \
      | "$BCFTOOLS" view --no-version --threads 2 \
            -T "^${LCR_BED}" \
            -Oz -o "${out}.tmp"
else
    echo "[$(date)] ${chrom}: filtering (exclude fails + indel<=50bp)"
    "$BCFTOOLS" view --no-version --threads 2 \
            -T "^${regions}" \
            -e 'abs(strlen(REF)-strlen(ALT))>50' \
            -Oz -o "${out}.tmp" "$in_bcf"
fi

mv "${out}.tmp" "$out"
"$TABIX" -p vcf -f "$out"

n=$("$BCFTOOLS" index -n "$out")
echo "[$(date)] ${chrom}: ${n} sites kept"
: > "$done_flag"
