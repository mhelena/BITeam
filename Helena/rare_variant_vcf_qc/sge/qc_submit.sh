#!/usr/bin/env bash
# Submits the five QC stages to SGE, each holding on the one before.
# Usage: ./qc_submit.sh [--redo chr10,chr16]    (settings in config.sh or $PIPELINE_CONFIG)
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
source "${PIPELINE_CONFIG:-${HERE}/config.sh}"

REDO=""
case "${1:-}" in
    --redo) REDO="${2:?--redo needs a comma-separated list of chromosomes}" ;;
    "")     ;;
    *)      echo "usage: $0 [--redo chr10,chr16]" >&2; exit 1 ;;
esac

case "$MODE" in
    allanc) export SAMPLE_FILE="" APPLY_HWE=0 APPLY_MEDIAN_DP_GQ=1 ;;
    eur)    export SAMPLE_FILE="${OUTDIR}/eur_samples.txt" APPLY_HWE=1 APPLY_MEDIAN_DP_GQ=0 ;;
    *)      echo "MODE must be allanc or eur, got '${MODE}'" >&2; exit 1 ;;
esac

export SCRIPTS="$HERE"
export HELPERS="$(cd "${HERE}/../helpers" && pwd)"

if [[ -n "${PYTHON_ENV_SOURCE:-}" ]]; then source "$PYTHON_ENV_SOURCE"; fi

PART="${OUTDIR}/_chunks"
LOG="${OUTDIR}/logs"
mkdir -p "$PART" "$LOG"

if [[ "$MODE" == eur ]]; then
    awk '{print $2}' "$EUR_KEEP" > "$SAMPLE_FILE"
    echo "EUR samples: $(wc -l < "$SAMPLE_FILE" | tr -d ' ')"
fi

# Built once and then reused, so array task numbers keep pointing at the same
# chromosome on reruns. Edit it before submitting to run a subset.
export CHROM_LIST="${OUTDIR}/chrom_list.txt"
if [[ ! -s "$CHROM_LIST" ]]; then
    "$BCFTOOLS" view -h "$RAW_VCF" 2>/dev/null | grep '^##contig' \
      | sed -E 's/.*ID=([^,>]+).*/\1/' \
      | grep -E '^(chr)?([0-9]+|X|Y|M|MT)$' > "$CHROM_LIST"
fi
N_CHR=$(wc -l < "$CHROM_LIST" | tr -d ' ')
echo "mode ${MODE}, chromosomes: ${N_CHR}"

# Finished chromosomes are skipped through their .done flags, so a plain rerun
# picks up where a failed run stopped. --redo also wipes the listed chromosomes
# (e.g. after an NFS error left a truncated BCF behind).
if [[ -n "$REDO" ]]; then
    for c in ${REDO//,/ }; do
        grep -qx "$c" "$CHROM_LIST" || { echo "$c is not in $CHROM_LIST" >&2; exit 1; }
        rm -f "${PART}/${c}.${MODE}.bcf" "${PART}/${c}.${MODE}.bcf.csi" "${PART}/${c}.subset.done" \
              "${PART}/${c}.ab_miss_fail.tsv" "${PART}/${c}.median_fail.tsv" "${PART}/${c}.hwe_fail.tsv" \
              "${PART}/${c}.ab_miss_stats.tsv.gz" "${PART}/${c}.fails.done" \
              "${PART}/${c}.qc.vcf.gz" "${PART}/${c}.qc.vcf.gz.tbi" "${PART}/${c}.filter.done"
        echo "cleared ${c}"
    done
fi

J="qc_${MODE}"

job0=$(qsub -terse -V \
            -N "${J}_0" -t 1-${N_CHR} \
            -o "${LOG}/${J}_0.\$TASK_ID.out" \
            -e "${LOG}/${J}_0.\$TASK_ID.err" \
            "${SCRIPTS}/qc_0_subset_chrom.sh" | cut -d. -f1)
echo "submitted step 0 (subset): ${job0}"

job1=$(qsub -terse -V \
            -N "${J}_1" -t 1-${N_CHR} \
            -hold_jid "${job0}" \
            -o "${LOG}/${J}_1.\$TASK_ID.out" \
            -e "${LOG}/${J}_1.\$TASK_ID.err" \
            "${SCRIPTS}/qc_1_fails_chrom.sh" | cut -d. -f1)
echo "submitted step 1 (fails) : ${job1}"

job2=$(qsub -terse -V \
            -N "${J}_2" \
            -hold_jid "${job1}" \
            -o "${LOG}/${J}_2.out" \
            -e "${LOG}/${J}_2.err" \
            "${SCRIPTS}/qc_2_union.sh" | cut -d. -f1)
echo "submitted step 2 (union) : ${job2}"

job3=$(qsub -terse -V \
            -N "${J}_3" -t 1-${N_CHR} \
            -hold_jid "${job2}" \
            -o "${LOG}/${J}_3.\$TASK_ID.out" \
            -e "${LOG}/${J}_3.\$TASK_ID.err" \
            "${SCRIPTS}/qc_3_filter_chrom.sh" | cut -d. -f1)
echo "submitted step 3 (filter): ${job3}"

job4=$(qsub -terse -V \
            -N "${J}_4" \
            -hold_jid "${job3}" \
            -o "${LOG}/${J}_4.out" \
            -e "${LOG}/${J}_4.err" \
            "${SCRIPTS}/qc_4_concat.sh" | cut -d. -f1)
echo "submitted step 4 (concat): ${job4}"

echo
echo "Final output will be: ${OUTDIR}/study_qc_${MODE}.vcf.gz"
echo "Fail-site audit:      ${OUTDIR}/study_qc_${MODE}.failsites.tsv"
