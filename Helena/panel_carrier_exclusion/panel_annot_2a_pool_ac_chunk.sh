#!/usr/bin/env bash
# Array job over aggV2 genomic_data chunks: pool AC, AN, AF, het and hom counts
# at the panel sites. Submitted by panel_annot_submit.sh.
#BSUB -q short
#BSUB -M 4000
#BSUB -n 2
#BSUB -R "rusage[mem=4000]"
set -euo pipefail
source "${PIPELINE_CONFIG:?set PIPELINE_CONFIG to the absolute path of config.sh}"
: "${GENO_CHUNKS_TXT:?}"

module load "$BCFTOOLS_MODULE"

WORK="${ANNOT_DIR}/_merged"
PART="${WORK}/_pool_ac"
mkdir -p "$PART"

regions="${WORK}/panel_annot_regions.tsv.gz"
[[ -s "$regions" && -s "${regions}.tbi" ]] || {
    echo "missing ${regions}(.tbi), run panel_annot_2pre_regions.sh first"; exit 1; }

vcf=$(sed -n "${LSB_JOBINDEX}p" "$GENO_CHUNKS_TXT")
[[ -s "$vcf" ]] || { echo "missing VCF at line ${LSB_JOBINDEX}: $vcf"; exit 1; }

base=$(basename "$vcf" .vcf.gz)
out="${PART}/${base}.tsv"
done_flag="${PART}/${base}.done"
empty_flag="${PART}/${base}.empty"

if [[ -e "$done_flag" || -e "$empty_flag" ]]; then
    echo "[$(date)] ${base}: already processed"; exit 0
fi

# fill-tags after -S so the counts are recomputed on the pool only
bcftools view --no-version --threads 2 \
        -R "$regions" \
        -S "$POOL_FILE" --force-samples \
        -Ou "$vcf" 2>/dev/null \
  | bcftools +fill-tags -Ou -- -t AC,AN,AF,AC_Het,AC_Hom \
  | bcftools query \
        -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AC\t%INFO/AN\t%INFO/AF\t%INFO/AC_Het\t%INFO/AC_Hom\n' \
        2>/dev/null > "$out" || true

n=$(wc -l < "$out" | tr -d ' ')
if [[ "$n" -eq 0 ]]; then
    rm -f "$out"
    : > "$empty_flag"
    echo "[$(date)] ${base}: no overlap with panel sites"
else
    : > "$done_flag"
    echo "[$(date)] ${base}: ${n} pool AC rows"
fi
