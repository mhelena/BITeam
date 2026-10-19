#!/bin/bash
# LDSC rg for every pair of traits, in trait table order.
set -euo pipefail
source "${PIPELINE_CONFIG:-$(dirname "$0")/config.sh}"

mkdir -p "$OUTDIR/rg"
names=($(cut -f1 "$OUTDIR/prep/traits.list"))  # macOS bash 3 has no mapfile

for ((i = 0; i < ${#names[@]}; i++)); do
  for ((j = i + 1; j < ${#names[@]}; j++)); do
    t1=${names[$i]}
    t2=${names[$j]}
    echo "rg $t1 vs $t2"
    $PYTHON "$LDSC_DIR/ldsc.py" \
      --rg "$OUTDIR/munged/$t1.sumstats.gz,$OUTDIR/munged/$t2.sumstats.gz" \
      --ref-ld-chr "$REF_LD." \
      --w-ld-chr "$W_LD." \
      --out "$OUTDIR/rg/rg_${t1}_${t2}" > /dev/null
  done
done
