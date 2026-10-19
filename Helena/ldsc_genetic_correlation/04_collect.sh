#!/bin/bash
# Pull rg, se, z, p and both observed-scale h2 out of the rg logs into one table.
set -euo pipefail
source "${PIPELINE_CONFIG:-$(dirname "$0")/config.sh}"

out="$OUTDIR/rg_summary.tsv"
names=($(cut -f1 "$OUTDIR/prep/traits.list"))  # macOS bash 3 has no mapfile

printf 'trait1\ttrait2\trg\tse\tz\tp\th2_obs1\th2_obs2\n' > "$out"
for ((i = 0; i < ${#names[@]}; i++)); do
  for ((j = i + 1; j < ${#names[@]}; j++)); do
    t1=${names[$i]}
    t2=${names[$j]}
    log="$OUTDIR/rg/rg_${t1}_${t2}.log"
    if [ ! -f "$log" ]; then
      echo "missing $log" >&2
      continue
    fi
    # columns after the summary header: p1 p2 rg se z p h2_obs h2_obs_se ...
    line=$(awk '/Summary of Genetic Correlation Results/ {f = 1; next}
                f && NF >= 12 && $1 != "p1" {print; exit}' "$log")
    if [ -z "$line" ]; then
      echo "no rg line in $log (check for errors)" >&2
      continue
    fi
    h1=$(grep -A2 "Heritability of phenotype 1" "$log" | awk '/Total Observed/ {print $5}')
    h2=$(grep -A2 "Heritability of phenotype 2" "$log" | awk '/Total Observed/ {print $5}')
    echo "$line" | awk -v a="$t1" -v b="$t2" -v h1="$h1" -v h2="$h2" \
      'BEGIN {OFS = "\t"} {print a, b, $3, $4, $5, $6, h1, h2}' >> "$out"
  done
done

column -t -s $'\t' "$out"
