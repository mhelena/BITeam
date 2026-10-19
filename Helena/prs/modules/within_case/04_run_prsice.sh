#!/usr/bin/env bash
# PRSice within the cases: every score in WITHIN_PRS_LIST against the binary outcome, and onset if present.
# Usage: bsub < 04_run_prsice.sh
#BSUB -P re_gecip_neurology
#BSUB -q long -M 16000 -n 2 -R "rusage[mem=16000]"
#BSUB -J within04 -o logs/%J_within04.out -e logs/%J_within04.err
set -euo pipefail
source "${PIPELINE_CONFIG:?export PIPELINE_CONFIG=/path/to/config.sh}"

TARGET="${WITHIN_DIR}/case_only.QC"
COV="${WITHIN_DIR}/cov.txt"

for f in "$PRSICE" "${TARGET}.bed" "${TARGET}.bim" "${TARGET}.fam" \
         "$COV" "${WITHIN_DIR}/pheno_binary.txt"; do
    [[ -e "$f" ]] || { echo "MISSING: $f" >&2; exit 1; }
done

run_prsice() {
    local pheno_file=$1 binary=$2 outdir=$3
    mkdir -p "$outdir"
    for phe in ${WITHIN_PRS_LIST}; do
        local var=GWAS_${phe}
        local base="${!var:-}"
        [[ -n "$base" && -e "$base" ]] || { echo "MISSING base for ${phe}: ${base}" >&2; exit 1; }
        echo "[$(date)] ${phe}, base ${base}, pheno $(basename ${pheno_file})"
        "$PRSICE" \
            --A1 A1 --A2 A2 \
            --all-score \
            --bar-levels 0.0001,0.001,0.01,0.05,0.1,0.5,1 \
            --base "$base" \
            --binary-target ${binary} \
            --bp BP --chr CHR \
            --fastscore \
            --ignore-fid \
            --model add \
            --pvalue P --se SE \
            --snp SNP --stat BETA \
            --cov-file  "$COV" \
            --pheno-file "$pheno_file" \
            --target "$TARGET" \
            --thread 2 \
            --print-snp \
            --out "${outdir}/PRS_${phe}"
    done
}

run_prsice "${WITHIN_DIR}/pheno_binary.txt" T "${WITHIN_DIR}/results_binary"

if [[ -s "${WITHIN_DIR}/pheno_onset.txt" ]]; then
    run_prsice "${WITHIN_DIR}/pheno_onset.txt" F "${WITHIN_DIR}/results_onset"
fi

ls -lh "${WITHIN_DIR}"/results_*/*.summary
