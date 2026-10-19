#!/usr/bin/env bash
# Optional: dump every panel-gene variant with its annotation and pool allele counts
# to ${ANNOT_DIR}/_merged/panel_annot_with_pool_ac.tsv.  Usage: ./panel_annot_submit.sh
set -euo pipefail

export PIPELINE_CONFIG="${PIPELINE_CONFIG:-$(cd "$(dirname "$0")" && pwd)/config.sh}"
source "$PIPELINE_CONFIG"

mkdir -p "$ANNOT_DIR" "$LOG_DIR"
[[ -s "$PANEL_GENES" ]] || { echo "missing $PANEL_GENES (panel_get_genes.R)"; exit 1; }
[[ -s "$POOL_FILE"   ]] || { echo "missing $POOL_FILE"; exit 1; }

export FA_CHUNKS_TXT="${ANNOT_DIR}/fa_chunks.txt"
ls "${FA_DIR}"/gel_mainProgramme_aggV2_*_annotated.vcf.gz > "$FA_CHUNKS_TXT"
n_fa=$(wc -l < "$FA_CHUNKS_TXT" | tr -d ' ')
echo "FA chunks: ${n_fa}"

export GENO_CHUNKS_TXT="${ANNOT_DIR}/geno_chunks.txt"
ls "${GENO_DIR}"/gel_mainProgramme_aggV2_*.vcf.gz > "$GENO_CHUNKS_TXT"
n_geno=$(wc -l < "$GENO_CHUNKS_TXT" | tr -d ' ')
echo "genomic_data chunks: ${n_geno}"

submit() {
    bsub -P "$LSF_PROJECT" -env "all" "$@" | awk '{print $2}' | tr -d '<>'
}

job1=$(submit -J "panel_annot1[1-${n_fa}]%200" \
       -o "${LOG_DIR}/annot1.%I.out" -e "${LOG_DIR}/annot1.%I.err" \
       < "${SCRIPT_DIR}/panel_annot_1_sites.sh")
echo "submitted annot 1 (sites): ${job1}"

job2pre=$(submit -J "panel_annot2pre" -w "done(${job1})" \
          -o "${LOG_DIR}/annot2pre.out" -e "${LOG_DIR}/annot2pre.err" \
          < "${SCRIPT_DIR}/panel_annot_2pre_regions.sh")
echo "submitted annot 2pre (regions): ${job2pre}"

job2a=$(submit -J "panel_annot2a[1-${n_geno}]%200" -w "done(${job2pre})" \
        -o "${LOG_DIR}/annot2a.%I.out" -e "${LOG_DIR}/annot2a.%I.err" \
        < "${SCRIPT_DIR}/panel_annot_2a_pool_ac_chunk.sh")
echo "submitted annot 2a (pool AC per chunk): ${job2a}"

job2b=$(submit -J "panel_annot2b" -w "done(${job2a})" \
        -o "${LOG_DIR}/annot2b.out" -e "${LOG_DIR}/annot2b.err" \
        < "${SCRIPT_DIR}/panel_annot_2b_pool_ac_concat.sh")
echo "submitted annot 2b (concat): ${job2b}"

job3=$(submit -J "panel_annot3" -w "done(${job2b})" \
       -o "${LOG_DIR}/annot3.out" -e "${LOG_DIR}/annot3.err" \
       < "${SCRIPT_DIR}/panel_annot_3_merge.sh")
echo "submitted annot 3 (merge): ${job3}"

echo "output: ${ANNOT_DIR}/_merged/panel_annot_with_pool_ac.tsv"
