#!/usr/bin/env bash
# Sensitivity run: longevity PRS again with a wide APOE window removed from the base GWAS.
# Usage: bsub < 03_prsice_longevity_no_apoe.sh
#BSUB -P re_gecip_neurology
#BSUB -q long -M 16000 -n 2 -R "rusage[mem=16000]"
#BSUB -J apoe03 -o logs/%J_apoe03.out -e logs/%J_apoe03.err
set -euo pipefail
source "${PIPELINE_CONFIG:?export PIPELINE_CONFIG=/path/to/config.sh}"

var=GWAS_${LONGEVITY_SCORE}
BASE="${!var:?no ${var} in config.sh}"
TARGET="${WORKDIR}/merged/case-HC-DEE-EPI.QC"
COV="${WORKDIR}/cov.txt"
PHENO="${WORKDIR}/pheno.txt"

OUT_DIR="${WORKDIR}/apoe"
BASE_NOAPOE="${OUT_DIR}/${LONGEVITY_SCORE}.noAPOE.tsv.gz"
OUT_PREFIX="${OUT_DIR}/PRS_${LONGEVITY_SCORE}_noAPOE"

mkdir -p "$OUT_DIR"

# hg38, about 500 kb either side of rs429358 (19:44908684), so the whole APOE LD
# block (TOMM40, APOC1) that co-tags the haplotypes goes too
APOE_CHR=19
APOE_START=44400000
APOE_END=45500000

echo "[$(date)] removing ${APOE_CHR}:${APOE_START}-${APOE_END} from ${BASE}"

zcat -f "$BASE" | awk -v FS='\t' -v OFS='\t' \
    -v chr="$APOE_CHR" -v start="$APOE_START" -v end="$APOE_END" '
    NR==1 {
        for (i=1; i<=NF; i++) col[$i]=i
        if (!("CHR" in col) || !("BP" in col)) {
            print "ERROR: CHR / BP columns missing in base sumstats" > "/dev/stderr"
            exit 1
        }
        print; next
    }
    {
        c = $col["CHR"]; sub(/^chr/, "", c)
        if (c == chr && $col["BP"] >= start && $col["BP"] <= end) next
        print
    }' | gzip > "$BASE_NOAPOE"

n_full=$(zcat -f "$BASE"     | tail -n +2 | wc -l)
n_kept=$(zcat "$BASE_NOAPOE" | tail -n +2 | wc -l)
echo "  SNPs in base:      ${n_full}"
echo "  SNPs after filter: ${n_kept}"
echo "  SNPs removed:      $(( n_full - n_kept ))"

"$PRSICE" \
    --A1 A1 --A2 A2 \
    --all-score \
    --bar-levels 0.0001,0.001,0.01,0.05,0.1,0.5,1 \
    --base "$BASE_NOAPOE" \
    --binary-target T \
    --bp BP --chr CHR \
    --fastscore \
    --ignore-fid \
    --model add \
    --pvalue P --se SE \
    --snp SNP --stat BETA \
    --cov-file  "$COV" \
    --pheno-file "$PHENO" \
    --target "$TARGET" \
    --thread 2 \
    --print-snp \
    --out "$OUT_PREFIX"

echo "[$(date)] done: ${OUT_PREFIX}.*"
