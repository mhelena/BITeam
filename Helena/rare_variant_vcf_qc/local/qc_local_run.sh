#!/usr/bin/env bash
# All-ancestry variant QC on one machine, in two passes over the whole VCF.
# Usage: ./qc_local_run.sh [--pass2-only]    (--pass2-only reuses pass 1's fail regions)
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
source "${PIPELINE_CONFIG:-${HERE}/config.sh}"
HELPER="${HERE}/../helpers/ab_miss_complete.py"

PASS2_ONLY=0
case "${1:-}" in
    --pass2-only) PASS2_ONLY=1 ;;
    "")           ;;
    *)            echo "usage: $0 [--pass2-only]" >&2; exit 1 ;;
esac

OUT="${OUTDIR}/${OUT_PREFIX}.vcf.gz"
LOG="${OUTDIR}/${OUT_PREFIX}.filters.log"
FAIL_TSV="${OUTDIR}/${OUT_PREFIX}.failsites.tsv"
STATS="${OUTDIR}/ab_miss_stats.tsv"
AB_FAIL="${OUTDIR}/ab_fail_regions.tsv"
FAIL_REGIONS="${OUTDIR}/all_fail_regions.tsv.gz"

mkdir -p "$OUTDIR"

if [[ "$PASS2_ONLY" == 1 ]]; then
    for f in "$IN" "$FAIL_REGIONS" "${FAIL_REGIONS}.tbi"; do
        [[ -e "$f" ]] || { echo "MISSING: $f" >&2; exit 1; }
    done
else
    for f in "$IN" "$HELPER"; do
        [[ -e "$f" ]] || { echo "MISSING: $f" >&2; exit 1; }
    done
    : > "$LOG"
fi

{
if [[ "$PASS2_ONLY" == 1 ]]; then echo "[$(date)] START (pass 2 only)"; else echo "[$(date)] START"; fi
echo "  IN     : $IN ($(du -h "$IN" | cut -f1))"
echo "  OUTDIR : $OUTDIR"

if [[ "$PASS2_ONLY" == 0 ]]; then
    N_SAMP=$("$BCFTOOLS" view -h "$IN" | tail -1 | awk '{print NF-9}')
    echo "  samples: $N_SAMP"

    echo "[$(date)] PASS 1: AB + miss + complete (the long one)"
    "$BCFTOOLS" query -f '%CHROM\t%POS\t%ID[\t%GT|%DP|%AD]\n' "$IN" \
      | "$PYTHON3" "$HELPER" "$STATS" "${OUTDIR}/ab_fail_ids.txt"

    awk -F'\t' 'NR>1 && $11!="PASS" {print $1"\t"$2}' "$STATS" \
      | sort -u -k1,1 -k2,2n > "$AB_FAIL"
    N_AB=$(wc -l < "$AB_FAIL" | tr -d ' ')
    echo "  AB/miss/complete fails: $N_AB sites"

    gzip -f "$STATS"
    rm -f "${OUTDIR}/ab_fail_ids.txt"

    # Only the AB/miss/complete fails go in the regions file; MEDIAN_DP/GQ and
    # indel length are applied as expressions in pass 2.
    sort -u -k1,1 -k2,2n "$AB_FAIL" \
      | "$BGZIP" -f > "$FAIL_REGIONS"
    "$TABIX" -s1 -b2 -e2 -f "$FAIL_REGIONS"

    awk '{print $1"\t"$2"\tAB_or_complete"}' "$AB_FAIL" \
      | sort -u -k1,1 -k2,2n > "$FAIL_TSV"
else
    echo "  FAIL_REGIONS: $FAIL_REGIONS ($(gzip -dc "$FAIL_REGIONS" | wc -l | tr -d ' ') sites)"
fi

echo "[$(date)] PASS 2: bcftools view (exclude AB fails + MEDIAN_DP/GQ + indel <= 50 bp)"
"$BCFTOOLS" view --no-version \
    -T "^${FAIL_REGIONS}" \
    -e 'INFO/MEDIAN_DP<10 || INFO/MEDIAN_GQ<15 || abs(strlen(REF)-strlen(ALT))>50' \
    -Oz -o "$OUT" "$IN"

"$TABIX" -p vcf -f "$OUT"

N_OUT=$("$BCFTOOLS" index -n "$OUT")
echo "[$(date)] FINAL kept sites: $N_OUT"
echo "         output: $OUT (+ .tbi)"
echo "         audit:  $FAIL_TSV"
echo "         chromosome coverage:"
"$TABIX" -l "$OUT" | sed 's/^/           /'
echo "[$(date)] DONE"
} 2>&1 | tee -a "$LOG"
