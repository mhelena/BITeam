#!/usr/bin/env bash
#BSUB -M 4000
#BSUB -n 1
# Canonical protein-coding annotation (HGVS, REVEL, SIFT, PolyPhen2 HVAR) for each
# carrier variant. Submitted by submit.sh annotate.
set -euo pipefail
source "${PIPELINE_CONFIG:?export PIPELINE_CONFIG=/absolute/path/to/config.sh}"

module load bcftools/1.16
mkdir -p "$WORK"

OUT="$VAR_ANNOT"

awk 'BEGIN{OFS="\t"} NR>0 && $1 ~ /^chr/ {print $1, $2-1, $2}' "$VARIANTS" \
  | sort -u > "${WORK}/variants.bed"
n_var=$(wc -l < "${WORK}/variants.bed" | tr -d ' ')
echo "[$(date)] input unique positions: ${n_var}"

printf "CHROM\tPOS\tREF\tALT\tSYMBOL\tFeature\tCANONICAL\tMANE_SELECT\tBIOTYPE\tConsequence\tHGVSc\tHGVSp\tREVEL\tSIFT_pred\tSIFT_score\tPolyphen2_HVAR_pred\tPolyphen2_HVAR_score\n" \
    > "$OUT"

# Rows that are missense only on an NMD transcript are dropped by the filter below.
if [[ -n "$TRANSCRIPT" ]]; then
    tx_filter="Feature==\"${TRANSCRIPT}\""
else
    tx_filter='CANONICAL=="YES"'
fi
filter="SYMBOL==\"${GENE}\" && BIOTYPE==\"protein_coding\" && ${tx_filter} && Consequence!~\"NMD_transcript_variant\""

for f in "$FA_DIR"/gel_mainProgramme_aggV2_"${REGION_CHR}"_*_annotated.vcf.gz; do
    base=$(basename "$f" _annotated.vcf.gz)
    end=${base##*_}; rest=${base%_*}; start=${rest##*_}
    if [[ "$end" -lt "$REGION_START" || "$start" -gt "$REGION_END" ]]; then continue; fi

    bcftools +split-vep \
        -R "${WORK}/variants.bed" \
        -c Feature,SYMBOL,BIOTYPE,CANONICAL,MANE_SELECT,Consequence,HGVSc,HGVSp,REVEL,SIFT_pred,SIFT_score,Polyphen2_HVAR_pred,Polyphen2_HVAR_score \
        -d \
        -f '%CHROM\t%POS\t%REF\t%ALT\t%SYMBOL\t%Feature\t%CANONICAL\t%MANE_SELECT\t%BIOTYPE\t%Consequence\t%HGVSc\t%HGVSp\t%REVEL\t%SIFT_pred\t%SIFT_score\t%Polyphen2_HVAR_pred\t%Polyphen2_HVAR_score\n' \
        -i "$filter" \
        "$f" >> "$OUT"
done

# Overlapping FA chunks can repeat a variant.
{ head -n1 "$OUT"; tail -n +2 "$OUT" | sort -u -k1,1 -k2,2n -k3,3 -k4,4; } \
  > "${OUT}.tmp" && mv "${OUT}.tmp" "$OUT"

n_out=$(( $(wc -l < "$OUT") - 1 ))
echo "[$(date)] annotated variants: ${n_out} -> ${OUT}"
echo "Positions missing from the output had no canonical protein-coding hit (e.g. NMD-only) or are not in the FA data."
