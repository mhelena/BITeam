#!/usr/bin/env bash
# Merges per-sample VCFs into one multi-sample VCF (bgzips and indexes inputs if needed).
# Usage: bash 00_merge_vcfs.sh <vcf_dir | file_with_one_vcf_path_per_line> <out.vcf.gz> [threads]
set -euo pipefail

INPUT=$1
OUT_VCF=$2
THREADS=${3:-4}

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

VCF_FILES=()
if [[ -d "$INPUT" ]]; then
    for f in "$INPUT"/*.vcf.gz "$INPUT"/*.vcf; do
        [[ -f "$f" ]] && VCF_FILES+=("$f")
    done
elif [[ -f "$INPUT" ]]; then
    while IFS= read -r line || [[ -n "$line" ]]; do
        [[ -z "$line" ]] && continue
        [[ ! -f "$line" ]] && { echo "ERROR: VCF not found: $line" >&2; exit 1; }
        VCF_FILES+=("$line")
    done < "$INPUT"
else
    echo "ERROR: $INPUT is neither a directory nor a file" >&2
    exit 1
fi

if [[ ${#VCF_FILES[@]} -eq 0 ]]; then
    echo "ERROR: No VCF files found" >&2
    exit 1
fi
log "Found ${#VCF_FILES[@]} VCF files"

INDEXED_FILES=()
for f in "${VCF_FILES[@]}"; do
    if [[ "$f" == *.vcf ]] && [[ ! -f "${f}.gz" ]]; then
        log "bgzipping $f"
        bgzip --threads "$THREADS" -c "$f" > "${f}.gz"
        f="${f}.gz"
    elif [[ "$f" == *.vcf ]]; then
        f="${f}.gz"
    fi
    if [[ ! -f "${f}.tbi" ]]; then
        tabix -p vcf "$f"
    fi
    INDEXED_FILES+=("$f")
done

TMPLIST=$(mktemp)
printf '%s\n' "${INDEXED_FILES[@]}" > "$TMPLIST"

# --missing-to-ref: a site absent from one sample's VCF is set to 0/0, not ./.
mkdir -p "$(dirname "$OUT_VCF")"
bcftools merge \
    --file-list "$TMPLIST" \
    --missing-to-ref \
    --threads "$THREADS" \
    -O z -o "$OUT_VCF"

bcftools index --tbi --threads "$THREADS" "$OUT_VCF"
rm -f "$TMPLIST"

N_SAMPLES=$(bcftools query -l "$OUT_VCF" | wc -l)
N_VARS=$(bcftools view -H "$OUT_VCF" | wc -l)
log "Merged VCF: $N_SAMPLES samples, $N_VARS variants: $OUT_VCF"
