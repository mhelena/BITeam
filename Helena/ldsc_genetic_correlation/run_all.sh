#!/bin/bash
# Runs the five steps in order. Set PIPELINE_CONFIG to use a config other than ./config.sh
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
export PIPELINE_CONFIG=${PIPELINE_CONFIG:-$here/config.sh}
source "$PIPELINE_CONFIG"

$PYTHON "$here/01_prep_sumstats.py" "$TRAITS" "$GWAS_DIR" "$OUTDIR" "$REF_LD"
bash "$here/02_munge.sh"
bash "$here/03_rg_all_pairs.sh"
bash "$here/04_collect.sh"
Rscript "$here/05_plot_rg.R" "$OUTDIR/rg_summary.tsv" "$TRAITS" "$OUTDIR" "$PLOT_FOCUS" "$PLOT_RG_LIMIT"
