#!/usr/bin/env bash
#$ -S /bin/bash
#$ -cwd
#$ -j n
#$ -l h_rt=02:00:00
#$ -l h_vmem=4G
#$ -l tmem=4G
# Merges the per-chromosome fail lists into one tabix-indexed exclusion file,
# plus an audit table of CHROM, POS and reason.
set -euo pipefail

: "${OUTDIR:?}"
: "${MODE:?}"
: "${APPLY_HWE:?}"
: "${APPLY_MEDIAN_DP_GQ:?}"

PART="${OUTDIR}/_chunks"
audit="${OUTDIR}/study_qc_${MODE}.failsites.tsv"
regions="${OUTDIR}/all_fail_regions.tsv.gz"

kinds="ab_miss_fail"
if [[ "$APPLY_MEDIAN_DP_GQ" == 1 ]]; then kinds="$kinds median_fail"; fi
if [[ "$APPLY_HWE" == 1 ]]; then kinds="$kinds hwe_fail"; fi

reason() {
    case "$1" in
        ab_miss_fail) echo AB_or_complete ;;
        median_fail)  echo MEDIAN_DP_GQ ;;
        hwe_fail)     echo HWE_EUR ;;
    esac
}

# sort -u on the key columns keeps one row per site, the first reason listed.
{
  for kind in $kinds; do
    r=$(reason "$kind")
    for f in "${PART}"/*."${kind}".tsv; do
      if [[ -s "$f" ]]; then awk -v r="$r" '{print $1"\t"$2"\t"r}' "$f"; fi
    done
  done
} | sort -u -k1,1 -k2,2n > "$audit"
n_audit=$(wc -l < "$audit" | tr -d ' ')
echo "[$(date)] audit rows: ${n_audit}"

for kind in $kinds; do
    cat "${PART}"/*."${kind}".tsv
done 2>/dev/null \
  | sort -u -k1,1 -k2,2n \
  | "$BGZIP" -f > "$regions"
"$TABIX" -s1 -b2 -e2 -f "$regions"

n_all=$(zcat "$regions" | wc -l | tr -d ' ')
echo "[$(date)] union fail sites: ${n_all}"
echo "         regions: $regions"
