#!/usr/bin/env bash
#BSUB -M 2000
#BSUB -n 1
#BSUB -R "rusage[mem=2000]"
# Lists missense sites in GENE that are absent from gnomAD, from the aggV2 VEP
# annotation (sites only, no genotypes). Submitted by submit.sh scan.
set -euo pipefail
source "${PIPELINE_CONFIG:?export PIPELINE_CONFIG=/absolute/path/to/config.sh}"

mkdir -p "$WORK"
module load bcftools/1.16

SITES="${WORK}/sites.tsv"

# FA chunk names carry their range: gel_mainProgramme_aggV2_<chr>_<start>_<end>_annotated.vcf.gz
fa_list="${WORK}/fa_inputs.txt"
: > "$fa_list"
for f in "${FA_DIR}"/gel_mainProgramme_aggV2_"${REGION_CHR}"_*_annotated.vcf.gz; do
    base=$(basename "$f" _annotated.vcf.gz)
    end=${base##*_}
    rest=${base%_*}
    start=${rest##*_}
    if [[ "$end" -ge "$REGION_START" && "$start" -le "$REGION_END" ]]; then
        echo "$f" >> "$fa_list"
    fi
done

n_fa=$(wc -l < "$fa_list" | tr -d ' ')
echo "[$(date)] FA chunks covering ${GENE}: ${n_fa}"
[[ "$n_fa" -gt 0 ]] || { echo "no FA chunk overlaps ${REGION}"; exit 1; }

# gnomADg_AF here is the gnomAD release bundled with the GEL annotation (v3.1.2 for
# VEP_109), which can lag the current gnomAD, so "absent" is relative to that release.
# bcftools 1.16 split-vep cannot compare a Float tag to "." in -i, hence the awk.
: > "$SITES"
while read -r fa; do
    bcftools +split-vep \
            -r "$REGION" \
            -c SYMBOL,Consequence,IMPACT,gnomADg_AF \
            -d \
            -f '%CHROM\t%POS\t%REF\t%ALT\t%SYMBOL\t%Consequence\t%IMPACT\t%gnomADg_AF\n' \
            -i "SYMBOL==\"${GENE}\" && Consequence~\"missense_variant\"" \
            "$fa" \
      | awk -F'\t' 'BEGIN{OFS="\t"} $8=="."' \
      | sort -u >> "$SITES"
done < "$fa_list"

sort -u -k1,1 -k2,2n -k3,3 -k4,4 "$SITES" -o "$SITES"

n_sites=$(wc -l < "$SITES" | tr -d ' ')
echo "[$(date)] ${GENE} missense sites absent from gnomAD: ${n_sites} -> ${SITES}"
