#!/usr/bin/env bash
#$ -S /bin/bash
#$ -cwd
#$ -j n
#$ -l h_rt=12:00:00
#$ -l h_vmem=6G
#$ -l tmem=6G
# Lists failing sites for one chromosome: allele balance, DP=0 missingness and
# completeness always; MEDIAN_DP/MEDIAN_GQ and HWE when switched on by MODE.
set -euo pipefail

: "${OUTDIR:?}"
: "${CHROM_LIST:?}"
: "${HELPERS:?}"
: "${MODE:?}"
: "${APPLY_HWE:?}"
: "${APPLY_MEDIAN_DP_GQ:?}"

PART="${OUTDIR}/_chunks"
chrom=$(sed -n "${SGE_TASK_ID}p" "$CHROM_LIST")
[[ -n "$chrom" ]] || { echo "no chrom at task ${SGE_TASK_ID}"; exit 1; }

in_bcf="${PART}/${chrom}.${MODE}.bcf"
[[ -s "$in_bcf" ]] || { echo "missing input: $in_bcf"; exit 1; }

ab_fail="${PART}/${chrom}.ab_miss_fail.tsv"
median_fail="${PART}/${chrom}.median_fail.tsv"
hwe_fail="${PART}/${chrom}.hwe_fail.tsv"
stats="${PART}/${chrom}.ab_miss_stats.tsv"
done_flag="${PART}/${chrom}.fails.done"

if [[ -e "$done_flag" ]]; then
    echo "[$(date)] ${chrom}: fails already computed"; exit 0
fi

echo "[$(date)] ${chrom}: AB + miss + complete"
"$BCFTOOLS" query -f '%CHROM\t%POS\t%ID[\t%GT|%DP|%AD]\n' "$in_bcf" \
  | "$PYTHON3" "${HELPERS}/ab_miss_complete.py" "$stats" "${PART}/${chrom}.ab_miss_fail_ids.txt"

awk -F'\t' 'NR>1 && $11!="PASS" {print $1"\t"$2}' "$stats" \
  | sort -u -k1,1 -k2,2n > "$ab_fail"
n_ab=$(wc -l < "$ab_fail" | tr -d ' ')
echo "  AB/miss/complete fails: ${n_ab}"

gzip -f "$stats"
rm -f "${PART}/${chrom}.ab_miss_fail_ids.txt"

if [[ "$APPLY_MEDIAN_DP_GQ" == 1 ]]; then
    # MEDIAN_DP and MEDIAN_GQ must already be in INFO. A missing value counts as a fail.
    echo "[$(date)] ${chrom}: MEDIAN_DP / MEDIAN_GQ"
    "$BCFTOOLS" query -f '%CHROM\t%POS\t%INFO/MEDIAN_DP\t%INFO/MEDIAN_GQ\n' "$in_bcf" \
      | awk -F'\t' '{
            dp = $3; gq = $4
            bad = 0
            if (dp == "." || dp+0 < 10) bad = 1
            if (gq == "." || gq+0 < 15) bad = 1
            if (bad) print $1"\t"$2
        }' \
      | sort -u -k1,1 -k2,2n > "$median_fail"
    n_med=$(wc -l < "$median_fail" | tr -d ' ')
    echo "  MEDIAN_DP<10 or MEDIAN_GQ<15 fails: ${n_med}"
fi

if [[ "$APPLY_HWE" == 1 ]]; then
    # bcftools +fill-tags only gained HWE in 1.16, so it is computed from GTs.
    echo "[$(date)] ${chrom}: HWE"
    "$BCFTOOLS" query -f '%CHROM\t%POS[\t%GT]\n' "$in_bcf" \
      | "$PYTHON3" "${HELPERS}/hwe_exact.py" 1e-5 \
      | sort -u -k1,1 -k2,2n > "$hwe_fail"
    n_hwe=$(wc -l < "$hwe_fail" | tr -d ' ')
    echo "  HWE<1e-5 fails: ${n_hwe}"
fi

: > "$done_flag"
echo "[$(date)] ${chrom}: done"
