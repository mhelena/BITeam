#!/usr/bin/env bash
# Build a merged hg38 CDS BED (PANEL_BED) for the panel genes from the GENCODE GTF.
# Usage: PIPELINE_CONFIG=/abs/path/config.sh bash panel_make_bed.sh
set -euo pipefail
source "${PIPELINE_CONFIG:?set PIPELINE_CONFIG to the absolute path of config.sh}"

GTF="$GENCODE_GTF"
OUT="$PANEL_BED"

[[ -s "$GTF" ]] || { echo "missing GTF $GTF, set GENCODE_GTF"; exit 1; }

command -v bedtools >/dev/null 2>&1 || module load "$BEDTOOLS_MODULE"

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

# Ensembl gene ids are preferred when available (panel_ensembl.txt, one ENSG
# per line) because they are stable across GENCODE releases; HGNC symbols drift.
if [[ -s "${WORK_DIR}/panel_ensembl.txt" ]]; then
    key="gene_id"
    awk 'NF>0 {sub(/\..*/,"",$1); print "gene_id \""$1"\";"}' \
        "${WORK_DIR}/panel_ensembl.txt" | sort -u > "$tmp/needles.txt"
    echo "[$(date)] keying on gene_id ($(wc -l < "$tmp/needles.txt") ENSGs)"
elif [[ -s "$PANEL_GENES" ]]; then
    key="gene_name"
    awk 'NF>0 {print "gene_name \""$1"\";"}' "$PANEL_GENES" \
      | sort -u > "$tmp/needles.txt"
    echo "[$(date)] keying on gene_name ($(wc -l < "$tmp/needles.txt") symbols)"
else
    echo "no panel_ensembl.txt or $PANEL_GENES"; exit 1
fi

zcat "$GTF" \
  | awk -F'\t' -v OFS='\t' -v needles="$tmp/needles.txt" -v key="$key" '
      BEGIN { while ((getline l < needles) > 0) names[l]=1 }
      $3 == "CDS" {
        info = $9
        if (key == "gene_id") {
            if (match(info, /gene_id "[^."]+(\.[0-9]+)?";/)) {
                tag = substr(info, RSTART, RLENGTH)
                # GENCODE gene_ids are versioned (ENSG...13), the needles are not
                sub(/\.[0-9]+";/, "\";", tag)
            } else next
        } else {
            if (match(info, /gene_name "[^"]+";/))
                tag = substr(info, RSTART, RLENGTH)
            else next
        }
        if (tag in names) print $1, $4-1, $5
      }' \
  | sort -k1,1 -k2,2n \
  | bedtools merge -i - > "$OUT"

n_bed=$(wc -l < "$OUT" | tr -d ' ')
echo "[$(date)] wrote ${n_bed} merged CDS intervals to $OUT"
