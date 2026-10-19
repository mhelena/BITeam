#!/usr/bin/env bash
#BSUB -M 3000
#BSUB -n 1
#BSUB -R "rusage[mem=3000]"
# Finds pool members carrying any of the sites from 02_sites.sh and writes the carrier
# list plus a carrier-by-variant table. Submitted by submit.sh scan.
set -euo pipefail
source "${PIPELINE_CONFIG:?export PIPELINE_CONFIG=/absolute/path/to/config.sh}"

module load bcftools/1.16
module load htslib/1.21

SITES="${WORK}/sites.tsv"
[[ -s "$SITES" ]] || { echo "no sites file from 02_sites.sh: $SITES"; exit 1; }
[[ -s "$POOL_FILE" ]] || { echo "missing pool file: $POOL_FILE"; exit 1; }

regions_tsv="${WORK}/regions.tsv.gz"
cut -f1,2 "$SITES" | sort -k1,1 -k2,2n -u | bgzip -f > "$regions_tsv"
tabix -s1 -b2 -e2 -f "$regions_tsv"

n_sites=$(zcat "$regions_tsv" | wc -l | tr -d ' ')
echo "[$(date)] querying ${n_sites} sites for ${POOL_FILE##*/}"

carriers="${WORK}/carriers_pool.txt"
raw_cv="${WORK}/carrier_variants_raw.tsv"
: > "$carriers"
: > "$raw_cv"

# Every chunk is tried: bcftools view -R returns nothing quickly when there is no overlap.
shopt -s nullglob
for vcf in "${GENO_DIR}"/gel_mainProgramme_aggV2_*.vcf.gz; do
    bcftools view --no-version \
            -R "$regions_tsv" \
            -S "$POOL_FILE" --force-samples \
            -Ou "$vcf" 2>/dev/null \
      | bcftools query \
            -f '[%SAMPLE\t%GT\t%CHROM\t%POS\t%REF\t%ALT\n]' \
            -i 'GT[*]="alt"' \
      | awk -F'\t' -v c="$carriers" -v t="$raw_cv" \
            '$2 ~ /1/ {print $1 >> c; print $1"\t"$3"\t"$4"\t"$5"\t"$6 >> t}' || true
done

sort -u -o "$carriers" "$carriers"
n_car=$(wc -l < "$carriers" | tr -d ' ')
echo "[$(date)] carriers in pool: ${n_car}"

sort -u "$POOL_FILE" > "${WORK}/pool_sorted.txt"
comm -12 "$carriers" "${WORK}/pool_sorted.txt" > "$CARRIERS"

# Regions are matched on position only, so keep only the exact site alleles here.
{
    printf "IID\tCHROM\tPOS\tREF\tALT\n"
    awk -F'\t' 'NR==FNR {s[$1 FS $2 FS $3 FS $4]; next}
                ($2 FS $3 FS $4 FS $5) in s' "$SITES" "$raw_cv" \
      | sort -u -k2,2 -k3,3n -k4,4 -k5,5 -k1,1
} > "$CARRIER_VARIANTS"
tail -n +2 "$CARRIER_VARIANTS" | cut -f2-5 | sort -u -k1,1 -k2,2n -k3,3 -k4,4 > "$VARIANTS"

n_final=$(wc -l < "$CARRIERS" | tr -d ' ')
echo "[$(date)] ${GENE} carriers: ${n_final} -> ${CARRIERS}"
echo "         carrier-by-variant table: ${CARRIER_VARIANTS}"
