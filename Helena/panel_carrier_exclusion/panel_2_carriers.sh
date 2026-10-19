#!/usr/bin/env bash
# Collect the per-chunk sites, find pool samples carrying any of them in aggV2
# genomic_data, and write FINAL_OUT = POOL_FILE minus carriers. Submitted by panel_submit.sh.
#BSUB -q short
#BSUB -M 3000
#BSUB -n 1
#BSUB -R "rusage[mem=3000]"
set -euo pipefail
source "${PIPELINE_CONFIG:?set PIPELINE_CONFIG to the absolute path of config.sh}"

module load "$BCFTOOLS_MODULE"
module load "$HTSLIB_MODULE"

mkdir -p "$FINAL_DIR"

all_sites="${WORK_DIR}/panel_sites_all.tsv"
cat "${SITES_DIR}"/*_sites.tsv 2>/dev/null \
  | sort -u -k1,1 -k2,2n -k3,3 -k4,4 > "$all_sites"

n_sites=$(wc -l < "$all_sites" | tr -d ' ')
echo "[$(date)] candidate panel sites: ${n_sites}"
[[ "$n_sites" -gt 0 ]] || { echo "no sites from panel_1_sites.sh"; exit 1; }

regions="${WORK_DIR}/panel_regions.tsv.gz"
cut -f1,2 "$all_sites" | sort -k1,1 -k2,2n -u | bgzip -f > "$regions"
tabix -s1 -b2 -e2 -f "$regions"

n_pool=$(wc -l < "$POOL_FILE" | tr -d ' ')
echo "[$(date)] querying ${n_sites} sites for ${n_pool} pool samples"

# Only query the genomic_data chunks that contain at least one site. Chunk
# files are named gel_mainProgramme_aggV2_<chrom>_<start>_<end>.vcf.gz.
chunks_to_query="${WORK_DIR}/chunks_to_query.txt"
: > "$chunks_to_query"

shopt -s nullglob
for chrom in $(zcat "$regions" | cut -f1 | sort -u); do
    read -r min_pos max_pos < <(zcat "$regions" \
        | awk -v c="$chrom" '$1==c { if (min=="" || $2<min) min=$2; if ($2>max) max=$2 }
                             END { print min, max }')

    for f in "$GENO_DIR"/gel_mainProgramme_aggV2_${chrom}_*.vcf.gz; do
        base=$(basename "$f" .vcf.gz)
        end=${base##*_}
        rest=${base%_*}
        start=${rest##*_}
        if [[ "$end" -ge "$min_pos" && "$start" -le "$max_pos" ]]; then
            n_in=$(zcat "$regions" | awk -v c="$chrom" -v s="$start" -v e="$end" \
                       '$1==c && $2>=s && $2<=e' | wc -l)
            if [[ "$n_in" -gt 0 ]]; then
                echo "$f" >> "$chunks_to_query"
            fi
        fi
    done
done
shopt -u nullglob

n_chunks=$(wc -l < "$chunks_to_query" | tr -d ' ')
echo "[$(date)] ${n_chunks} genomic_data chunks to query"

carriers="${WORK_DIR}/panel_carriers.txt"
: > "$carriers"

i=0
while read -r vcf; do
    i=$((i+1))
    echo "[$(date)]   chunk ${i}/${n_chunks}: $(basename "$vcf")"
    bcftools view --no-version \
            -R "$regions" \
            -S "$POOL_FILE" --force-samples \
            -Ou "$vcf" 2>/dev/null \
      | bcftools query \
            -f '[%SAMPLE\t%GT\n]' \
            -i 'GT[*]="alt"' 2>/dev/null \
      | awk '$2 ~ /1/ {print $1}' >> "$carriers" || true
done < "$chunks_to_query"

sort -u -o "$carriers" "$carriers"
n_car=$(wc -l < "$carriers" | tr -d ' ')
echo "[$(date)] pool carriers of panel sites absent from gnomAD: ${n_car}"
echo "         written to: $carriers"

sort -u "$POOL_FILE" > "${WORK_DIR}/pool_sorted.txt"
comm -23 "${WORK_DIR}/pool_sorted.txt" "$carriers" > "$FINAL_OUT"

n_final=$(wc -l < "$FINAL_OUT" | tr -d ' ')
echo "[$(date)] final controls (pool minus carriers): ${n_final}"
echo "         written to: $FINAL_OUT"
