#!/usr/bin/env bash
# Array job over aggV2 VEP chunks: dump every panel-gene variant row, with no
# gnomAD or impact filter. Submitted by panel_annot_submit.sh.
#BSUB -q short
#BSUB -M 2000
#BSUB -n 1
#BSUB -R "rusage[mem=2000]"
set -euo pipefail
source "${PIPELINE_CONFIG:?set PIPELINE_CONFIG to the absolute path of config.sh}"
: "${FA_CHUNKS_TXT:?}"

mkdir -p "$ANNOT_DIR"
module load "$BCFTOOLS_MODULE"

fa=$(sed -n "${LSB_JOBINDEX}p" "$FA_CHUNKS_TXT")
[[ -s "$fa" ]] || { echo "missing FA file at line ${LSB_JOBINDEX}: $fa"; exit 1; }

base=$(basename "$fa" .vcf.gz)
out="${ANNOT_DIR}/${base}_annot.tsv"
done_flag="${ANNOT_DIR}/${base}.annot.done"
empty_flag="${ANNOT_DIR}/${base}.annot.empty"

if [[ -e "$done_flag" || -e "$empty_flag" ]]; then
    echo "[$(date)] ${base}: already processed"; exit 0
fi

bcftools +split-vep \
        -c SYMBOL,Consequence,IMPACT,gnomADg_AF \
        -d \
        -f '%CHROM\t%POS\t%REF\t%ALT\t%SYMBOL\t%Consequence\t%IMPACT\t%gnomADg_AF\n' \
        "$fa" \
  | awk -F'\t' -v genes="$PANEL_GENES" '
        BEGIN { OFS="\t"; while ((getline g < genes) > 0) keep[g] = 1 }
        keep[$5] { print }
    ' \
  | sort -u > "$out"

n=$(wc -l < "$out" | tr -d ' ')
if [[ "$n" -eq 0 ]]; then
    rm -f "$out"
    : > "$empty_flag"
    echo "[$(date)] ${base}: no panel-gene rows"
else
    : > "$done_flag"
    echo "[$(date)] ${base}: ${n} annotation rows"
fi
