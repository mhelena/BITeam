#!/usr/bin/env bash
# Concatenate the per-chunk pool AC files from 2a into _merged/panel_pool_ac.tsv.
# Submitted by panel_annot_submit.sh.
#BSUB -q short
#BSUB -M 4000
#BSUB -n 1
#BSUB -R "rusage[mem=4000]"
set -euo pipefail
source "${PIPELINE_CONFIG:?set PIPELINE_CONFIG to the absolute path of config.sh}"

WORK="${ANNOT_DIR}/_merged"
PART="${WORK}/_pool_ac"
out="${WORK}/panel_pool_ac.tsv"

n_parts=$(ls "${PART}"/*.tsv 2>/dev/null | wc -l | tr -d ' ')
[[ "$n_parts" -gt 0 ]] || { echo "no per-chunk outputs in ${PART}"; exit 1; }
echo "[$(date)] concatenating ${n_parts} per-chunk pool AC files"

# a site can sit in two overlapping chunks
sort -u -k1,1 -k2,2n -k3,3 -k4,4 "${PART}"/*.tsv > "$out"

n=$(wc -l < "$out" | tr -d ' ')
echo "[$(date)] pool AC rows: ${n}"
echo "         written to: $out"
