#!/usr/bin/env bash
# Submits the three extraction jobs for one cohort, each waiting for the previous one.
# Usage: ./submit.sh [queue]    (default queue: long)
set -euo pipefail
: "${PIPELINE_CONFIG:?export PIPELINE_CONFIG=/absolute/path/to/config.sh first}"
[[ "$PIPELINE_CONFIG" == /* ]] || { echo "PIPELINE_CONFIG must be an absolute path"; exit 1; }
source "$PIPELINE_CONFIG"

QUEUE="${1:-long}"
HERE="$(cd "$(dirname "$0")" && pwd)"

LOG=${WORKDIR}/logs
export SAMPLES=${WORKDIR}/samples_eur08.txt
export OUTDIR=${SCRATCH}/chunks
export CHRDIR=${SCRATCH}/perchr
export FINALDIR=${WORKDIR}/final
export AB_PY=${HERE}/_ab_miss_complete.py
mkdir -p "$LOG" "$OUTDIR" "$CHRDIR" "$FINALDIR"

[[ -s "$SAMPLES" ]]    || { echo "missing ${SAMPLES}, run 01_eur_filter.sh first"; exit 1; }
[[ -s "$CHUNKS_TXT" ]] || { echo "missing ${CHUNKS_TXT}"; exit 1; }
[[ -s "$AB_PY" ]]      || { echo "missing ${AB_PY}"; exit 1; }

N_CHUNKS=$(wc -l < "$CHUNKS_TXT" | tr -d ' ')
echo "[$(date)] ${COHORT}: $(wc -l < "$SAMPLES" | tr -d ' ') samples, ${N_CHUNKS} chunks, queue ${QUEUE}"

# Variables must be exported before bsub: the GEL esub fails on "-env all,VAR=val".
J1="subset_${COHORT}"
bsub -P "$LSF_PROJECT" -J "${J1}[1-${N_CHUNKS}]%500" -q "$QUEUE" \
     -o "${LOG}/s1.%J.%I.out" -e "${LOG}/s1.%J.%I.err" \
     -env "all" < "${HERE}/1_subset_plink.sh"

J2="perchr_${COHORT}"
bsub -P "$LSF_PROJECT" -J "${J2}[1-22]" -q "$QUEUE" -w "ended(${J1})" \
     -o "${LOG}/s2.%J.%I.out" -e "${LOG}/s2.%J.%I.err" \
     -env "all" < "${HERE}/2_per_chr_merge.sh"

bsub -P "$LSF_PROJECT" -J "merge_all_${COHORT}" -q "$QUEUE" -w "ended(${J2})" \
     -o "${LOG}/s3.%J.out" -e "${LOG}/s3.%J.err" \
     -env "all" < "${HERE}/3_all_chr_merge.sh"

echo "[$(date)] submitted. Logs in ${LOG}/"
echo "final BED will be ${FINALDIR}/${COHORT}_ALL.{bed,bim,fam}"
