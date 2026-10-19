#!/usr/bin/env bash
# Concatenate the per-chunk annotation rows and build the tabix-indexed regions
# file used by the 2a array. Submitted by panel_annot_submit.sh.
#BSUB -q short
#BSUB -M 4000
#BSUB -n 1
#BSUB -R "rusage[mem=4000]"
set -euo pipefail
source "${PIPELINE_CONFIG:?set PIPELINE_CONFIG to the absolute path of config.sh}"

module load "$HTSLIB_MODULE"

WORK="${ANNOT_DIR}/_merged"
mkdir -p "$WORK"

all_annot="${WORK}/panel_annot_all.tsv"
cat "${ANNOT_DIR}"/*_annot.tsv 2>/dev/null \
  | sort -u -k1,1 -k2,2n -k3,3 -k4,4 > "$all_annot"

n_rows=$(wc -l < "$all_annot" | tr -d ' ')
echo "[$(date)] annotation rows: ${n_rows}"
[[ "$n_rows" -gt 0 ]] || { echo "no annotation rows"; exit 1; }

regions="${WORK}/panel_annot_regions.tsv.gz"
cut -f1,2 "$all_annot" | sort -k1,1 -k2,2n -u | bgzip -f > "$regions"
tabix -s1 -b2 -e2 -f "$regions"

n_sites=$(zcat "$regions" | wc -l | tr -d ' ')
echo "[$(date)] unique sites to query in pool: ${n_sites}"
echo "         regions written to: $regions"
