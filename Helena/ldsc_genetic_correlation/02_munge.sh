#!/bin/bash
# Run munge_sumstats.py on every prepared trait. N flags come from prep/traits.list.
set -euo pipefail
source "${PIPELINE_CONFIG:-$(dirname "$0")/config.sh}"

mkdir -p "$OUTDIR/munged"
extra=""
[ -n "$MERGE_ALLELES" ] && extra="--merge-alleles $MERGE_ALLELES"

while IFS=$'\t' read -r name nargs; do
  echo "munging $name ($nargs)"
  # shellcheck disable=SC2086
  $PYTHON "$LDSC_DIR/munge_sumstats.py" \
    --sumstats "$OUTDIR/prep/$name.txt" \
    $nargs $extra \
    --snp SNP --a1 A1 --a2 A2 --p P --signed-sumstats BETA,0 \
    --out "$OUTDIR/munged/$name" > /dev/null
  grep -E "^Writing summary statistics|Mean chi\^2" "$OUTDIR/munged/$name.log" || true
done < "$OUTDIR/prep/traits.list"
