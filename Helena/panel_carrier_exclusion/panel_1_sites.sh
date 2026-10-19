#!/usr/bin/env bash
# Array job over aggV2 VEP chunks: write panel sites with HIGH/MODERATE impact
# that are absent from gnomADg v3.1.2. Submitted by panel_submit.sh.
#BSUB -q short
#BSUB -M 1500
#BSUB -n 1
#BSUB -R "rusage[mem=1500]"
set -euo pipefail
source "${PIPELINE_CONFIG:?set PIPELINE_CONFIG to the absolute path of config.sh}"
: "${FA_CHUNKS_TXT:?}"

mkdir -p "$SITES_DIR"
module load "$BCFTOOLS_MODULE"

fa=$(sed -n "${LSB_JOBINDEX}p" "$FA_CHUNKS_TXT")
[[ -s "$fa" ]] || { echo "missing FA file at line ${LSB_JOBINDEX}: $fa"; exit 1; }

base=$(basename "$fa" .vcf.gz)
out="${SITES_DIR}/${base}_sites.tsv"
done_flag="${SITES_DIR}/${base}.done"
empty_flag="${SITES_DIR}/${base}.empty"

if [[ -e "$done_flag" || -e "$empty_flag" ]]; then
    echo "[$(date)] ${base}: already processed"; exit 0
fi

region_opt=()
genes_file=""
if [[ "$PANEL_MODE" == bed ]]; then
    region_opt=(-R "$PANEL_BED")
else
    genes_file="$PANEL_GENES"
fi

# -d gives one row per transcript so panel-gene hits on non-canonical isoforms
# are not lost. gnomAD absence is tested in awk because split-vep -i in
# bcftools 1.16 will not compare the Float gnomADg_AF to ".".
bcftools +split-vep \
        ${region_opt[@]+"${region_opt[@]}"} \
        -c SYMBOL,Consequence,IMPACT,gnomADg_AF \
        -d \
        -f '%CHROM\t%POS\t%REF\t%ALT\t%SYMBOL\t%Consequence\t%IMPACT\t%gnomADg_AF\n' \
        "$fa" \
  | awk -F'\t' -v genes="$genes_file" '
        BEGIN { OFS="\t"; if (genes != "") while ((getline g < genes) > 0) keep[g] = 1 }
        (genes == "" || keep[$5]) && $8=="." && ($7=="HIGH" || $7=="MODERATE") { print }
    ' \
  | sort -u > "$out"

n=$(wc -l < "$out" | tr -d ' ')
if [[ "$n" -eq 0 ]]; then
    rm -f "$out"
    : > "$empty_flag"
    echo "[$(date)] ${base}: no panel sites absent from gnomAD"
else
    : > "$done_flag"
    echo "[$(date)] ${base}: ${n} candidate sites"
fi
