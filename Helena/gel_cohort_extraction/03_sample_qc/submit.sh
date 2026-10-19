#!/usr/bin/env bash
# Submits sample QC (Ti/Tv outliers, then HWE) for one cohort once 02_extract_plink has finished.
# Usage: ./submit.sh [queue]    (default queue: short)
set -euo pipefail
: "${PIPELINE_CONFIG:?export PIPELINE_CONFIG=/absolute/path/to/config.sh first}"
[[ "$PIPELINE_CONFIG" == /* ]] || { echo "PIPELINE_CONFIG must be an absolute path"; exit 1; }
source "$PIPELINE_CONFIG"

QUEUE="${1:-short}"
HERE="$(cd "$(dirname "$0")" && pwd)"
LOG=${WORKDIR}/logs
mkdir -p "$LOG"

IN=${WORKDIR}/final/${COHORT}_ALL
[[ -s "${IN}.bed" ]] || { echo "missing ${IN}.bed, run 02_extract_plink/submit.sh first"; exit 1; }

bsub -P "$LSF_PROJECT" -J "sampleqc_${COHORT}" -q "$QUEUE" \
     -o "${LOG}/sampleqc.%J.out" -e "${LOG}/sampleqc.%J.err" \
     -env "all" < "${HERE}/sample_qc.sh"

echo "QC'd BED will be ${WORKDIR}/qc/${COHORT}_qc.{bed,bim,fam}"
