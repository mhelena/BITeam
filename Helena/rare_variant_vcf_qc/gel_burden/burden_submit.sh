#!/usr/bin/env bash
# build: one VCF per cohort in COHORTS (in parallel), then a merge job that waits for all.
# acgt0: drop AC=0 sites from ACGT0_IN after variant QC.  Usage: ./burden_submit.sh build|acgt0
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
source "${PIPELINE_CONFIG:-${HERE}/config.sh}"

LOG="${OUTDIR}/logs"
mkdir -p "$LOG"

case "${1:-}" in
build)
    # Variables must be exported before bsub: the GEL esub fails on "-env all,VAR=val".
    JOBS=""
    for C in $COHORTS; do
        export COHORT="$C"
        export PERCHR="$(perchr_dir "$C")"
        if ! ls "${PERCHR}"/${C}_chr*.bed >/dev/null 2>&1; then
            echo "ERROR: no per-chr BEDs found at ${PERCHR}/${C}_chr*.bed"
            exit 1
        fi
        JOB=$(bsub -P "$LSF_PROJECT" \
            -J "burden_vcf_${C}" \
            -o "${LOG}/burden_vcf_${C}.out" \
            -e "${LOG}/burden_vcf_${C}.err" \
            -env "all" \
            < "${HERE}/burden_1_make_cohort_vcf.sh" \
            | awk '{gsub(/[<>]/,""); print $2}')
        echo "submitted ${C} -> jobid ${JOB}"
        JOBS="$JOBS $JOB"
    done

    HOLD=""
    for j in $JOBS; do HOLD="${HOLD:+${HOLD}&&}done(${j})"; done
    JOB2=$(bsub -P "$LSF_PROJECT" \
        -J "burden_vcf_merge_${MERGE_OUT_TAG}" \
        -w "${HOLD}" \
        -o "${LOG}/burden_vcf_merge_${MERGE_OUT_TAG}.out" \
        -e "${LOG}/burden_vcf_merge_${MERGE_OUT_TAG}.err" \
        -env "all" \
        < "${HERE}/burden_2_merge.sh" \
        | awk '{gsub(/[<>]/,""); print $2}')
    echo "submitted merge -> jobid ${JOB2}  (holds for${JOBS})"

    echo
    echo "Final output will be: ${OUTDIR}/${MERGE_OUT_TAG}_for_burden.vcf.gz"
    echo "Watch with:  bjobs -a -u \$USER | grep burden_vcf"
    ;;
acgt0)
    JOB_ID=$(bsub -q "$ACGT0_QUEUE" -P "$LSF_PROJECT" -n 2 -M 16000 \
                  -R "rusage[mem=16000] span[hosts=1]" \
                  -W 24:00 \
                  -o "${LOG}/ACgt0.%J.out" -e "${LOG}/ACgt0.%J.err" \
                  -J burden_ACgt0 \
                  bash "${HERE}/burden_3_acgt0.sh" |
             awk '{print $2}' | tr -d '<>')

    echo "Submitted job ${JOB_ID}"
    echo "Watch with:  bjobs ${JOB_ID}"
    echo "Tail log:    tail -f ${LOG}/ACgt0.${JOB_ID}.out"
    ;;
*)
    echo "usage: $0 build|acgt0" >&2
    exit 1
    ;;
esac
