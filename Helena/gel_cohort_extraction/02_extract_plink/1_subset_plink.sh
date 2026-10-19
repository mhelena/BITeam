#!/usr/bin/env bash
#BSUB -q medium
#BSUB -M 1500
#BSUB -n 1
#BSUB -R "rusage[mem=1500]"
# One aggV2 chunk per array index: site QC, cohort subset, AB/missingness/completeness
# filter, then a PLINK BED. Submitted by submit.sh, which sets SAMPLES, OUTDIR and AB_PY.
set -euo pipefail
: "${PIPELINE_CONFIG:?export PIPELINE_CONFIG=/absolute/path/to/config.sh first}"
source "$PIPELINE_CONFIG"

: "${COHORT:?}"
: "${SAMPLES:?}"
: "${OUTDIR:?}"
: "${CHUNKS_TXT:?}"
: "${AB_PY:?}"

mkdir -p "$OUTDIR"
export TMPDIR="$OUTDIR"

module load "$BCFTOOLS_MODULE"
module load "$PLINK2_MODULE"

i=$(sed -n "${LSB_JOBINDEX}p" "$CHUNKS_TXT" | awk '{print $1}')
in_vcf="${AGGV2_DIR}/gel_mainProgramme_aggV2_${i}.vcf.gz"
prefix="${OUTDIR}/chunk_${i}"

if [[ -s "${prefix}.bed" && -s "${prefix}.bim" && -s "${prefix}.fam" ]] \
   || [[ -e "${prefix}.empty" ]]; then
    echo "[$(date)] chunk ${i}: already done, skipping"; exit 0
fi

REGIONS=$(seq 1 22 | sed 's/^/chr/' | tr '\n' ',' | sed 's/,$//')
tmp_vcf="${OUTDIR}/.tmp_${i}.vcf.gz"
site_qc="${OUTDIR}/chunk_${i}_siteqc.tsv"
fail_lst="${OUTDIR}/chunk_${i}_fail.txt"
fail_bed="${OUTDIR}/chunk_${i}_fail.bed"

# aggV2 site metrics are camelCase INFO tags (medianDepthNonMiss, medianGQ).
bcftools view \
        --no-version --threads 1 \
        --targets "$REGIONS" \
        --types snps \
        --min-alleles 2 --max-alleles 2 \
        --apply-filters PASS \
        -S "$SAMPLES" --force-samples \
        "$in_vcf" -Ou \
  | bcftools filter \
        --no-version --threads 1 \
        --exclude 'INFO/medianDepthNonMiss<10 || INFO/medianGQ<15' -Ou \
  | bcftools +fill-tags --threads 1 -- -t AC,AN,AF \
  | bcftools +missing2ref --threads 1 \
  | bcftools view --threads 1 -i 'INFO/AC>0' -Oz -o "$tmp_vcf"

bcftools index --threads 1 "$tmp_vcf"

# Count from the index; piping into head would trip pipefail with SIGPIPE.
n_var=$(bcftools index --stats "$tmp_vcf" | awk '{s+=$3} END{print s+0}')
if [[ "$n_var" -eq 0 ]]; then
    echo "[$(date)] chunk ${i}: empty after site+cohort QC"
    rm -f "$tmp_vcf" "${tmp_vcf}.csi"
    : > "${prefix}.empty"
    exit 0
fi

bcftools query \
        --format '%CHROM\t%POS\t%ID[\t%GT|%DP|%AD]\n' \
        "$tmp_vcf" \
  | python3 "$AB_PY" "$site_qc" "$fail_lst"

awk 'BEGIN{OFS="\t"} NR>1 && $11!="PASS" {print $1, $2-1, $2}' \
    "$site_qc" > "$fail_bed"
n_fail=$(wc -l < "$fail_bed" | tr -d ' ')
echo "[$(date)] chunk ${i}: ${n_fail} sites failing AB/miss/comp"

EXCLUDE_OPT=""
if [[ "$n_fail" -gt 0 ]]; then
    EXCLUDE_OPT="--exclude bed1 ${fail_bed}"
fi

# MAF, geno and palindrome filters need the whole cohort, so they wait for 3_all_chr_merge.sh.
# shellcheck disable=SC2086
plink2 \
    --vcf "$tmp_vcf" \
    $EXCLUDE_OPT \
    --snps-only just-acgt \
    --max-alleles 2 \
    --autosome \
    --vcf-half-call m \
    --set-all-var-ids '@:#:$r:$a' \
    --new-id-max-allele-len 60 missing \
    --double-id \
    --real-ref-alleles \
    --memory 1200 \
    --threads 1 \
    --make-bed \
    --out "$prefix"

rm -f "$tmp_vcf" "${tmp_vcf}.csi"
echo "[$(date)] chunk ${i}: done"
