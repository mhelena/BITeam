#!/usr/bin/env bash
# Submit the carrier scan: per-chunk site discovery (array), then carrier lookup
# and FINAL_OUT = POOL_FILE minus carriers.  Usage: ./panel_submit.sh
set -euo pipefail

export PIPELINE_CONFIG="${PIPELINE_CONFIG:-$(cd "$(dirname "$0")" && pwd)/config.sh}"
source "$PIPELINE_CONFIG"

mkdir -p "$WORK_DIR" "$SITES_DIR" "$FINAL_DIR" "$LOG_DIR"

[[ -s "$POOL_FILE" ]] || { echo "missing $POOL_FILE"; exit 1; }
case "$PANEL_MODE" in
    bed)    [[ -s "$PANEL_BED"   ]] || { echo "missing $PANEL_BED (panel_make_bed.sh)"; exit 1; } ;;
    symbol) [[ -s "$PANEL_GENES" ]] || { echo "missing $PANEL_GENES (panel_get_genes.R)"; exit 1; } ;;
    *)      echo "PANEL_MODE must be bed or symbol, got '$PANEL_MODE'"; exit 1 ;;
esac

export FA_CHUNKS_TXT="${WORK_DIR}/fa_chunks.txt"
if [[ ! -s "$FA_CHUNKS_TXT" ]]; then
    ls "${FA_DIR}"/gel_mainProgramme_aggV2_*_annotated.vcf.gz \
      | sort -V > "$FA_CHUNKS_TXT"
fi
N_CHUNKS=$(wc -l < "$FA_CHUNKS_TXT" | tr -d ' ')
[[ "$N_CHUNKS" -gt 0 ]] || { echo "no FA chunks listed in $FA_CHUNKS_TXT"; exit 1; }

echo "[$(date)] panel scan (${PANEL_MODE}): ${N_CHUNKS} FA chunks, $(wc -l < "$POOL_FILE") candidates"

J1="panel_sites"
bsub -P "$LSF_PROJECT" -J "${J1}[1-${N_CHUNKS}]%500" \
     -o "${LOG_DIR}/panel_s1.%J.%I.out" -e "${LOG_DIR}/panel_s1.%J.%I.err" \
     -env "all" \
     < "${SCRIPT_DIR}/panel_1_sites.sh"

bsub -P "$LSF_PROJECT" -J "panel_carriers" \
     -w "ended(${J1})" \
     -o "${LOG_DIR}/panel_s2.%J.out" -e "${LOG_DIR}/panel_s2.%J.err" \
     -env "all" \
     < "${SCRIPT_DIR}/panel_2_carriers.sh"

echo "[$(date)] both stages submitted, final list will be written to $FINAL_OUT"
