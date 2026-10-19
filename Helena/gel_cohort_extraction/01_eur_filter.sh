#!/usr/bin/env bash
# Keeps samples with predicted European ancestry >= 0.8 and writes their platekeys.
# Needs ${WORKDIR}/samples.txt (platekeys or participant IDs, one per line). Run: bash 01_eur_filter.sh
set -euo pipefail
: "${PIPELINE_CONFIG:?export PIPELINE_CONFIG=/absolute/path/to/config.sh first}"
source "$PIPELINE_CONFIG"

cd "$WORKDIR"
[[ -s samples.txt ]] || { echo "Put your sample list in ${WORKDIR}/samples.txt first"; exit 1; }

awk -F'\t' 'NR==FNR {sub(/\r$/, ""); if ($1 != "") want[$1]; next}
            FNR > 1 && (($1 in want) || ($2 in want)) && $3+0 >= 0.8 {print $1}' \
    samples.txt "$ANCESTRY_SEX_TSV" | sort -u > samples_eur08.txt

n_in=$(awk 'NF' samples.txt | sort -u | wc -l | tr -d ' ')
n_found=$(awk -F'\t' 'NR==FNR {sub(/\r$/, ""); want[$1]; next}
                      FNR > 1 && (($1 in want) || ($2 in want))' \
          samples.txt "$ANCESTRY_SEX_TSV" | wc -l | tr -d ' ')
n_eur=$(wc -l < samples_eur08.txt | tr -d ' ')

echo "IDs in samples.txt          : ${n_in}"
echo "found in GEL ancestry table : ${n_found}"
echo "EUR >= 0.8, kept            : ${n_eur}   (${WORKDIR}/samples_eur08.txt, platekeys)"
