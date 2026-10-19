#!/usr/bin/env bash
# Submits the LSF steps: "scan" runs 02_sites.sh then 03_carriers.sh, "annotate" runs 06.
# Usage: export PIPELINE_CONFIG=/absolute/path/to/config.sh; bash submit.sh scan|annotate
set -euo pipefail
source "${PIPELINE_CONFIG:?export PIPELINE_CONFIG=/absolute/path/to/config.sh}"
here=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$OUT_DIR" "$WORK" "$LOG_DIR"

# bsub -env "all,VAR=val" breaks the RE esub parser, so rely on exported vars and -env all.
sub() {
    local name=$1 script=$2; shift 2
    bsub -P "$LSF_PROJECT" -q "$LSF_QUEUE" -J "$name" -env "all" \
         -o "${LOG_DIR}/${name}.%J.out" -e "${LOG_DIR}/${name}.%J.err" \
         "$@" < "${here}/${script}"
}

case "${1:-}" in
    scan)
        [[ -s "$POOL_FILE" ]] || { echo "missing $POOL_FILE, run 01_build_pool.R first"; exit 1; }
        sub "${GENE}_sites" 02_sites.sh
        sub "${GENE}_carriers" 03_carriers.sh -w "ended(${GENE}_sites)"
        echo "carrier list will be written to ${CARRIERS}"
        ;;
    annotate)
        [[ -s "$VARIANTS" ]] || { echo "missing $VARIANTS, run the scan first"; exit 1; }
        sub "${GENE}_annotate" 06_annotate_variants.sh
        ;;
    *)
        echo "usage: bash submit.sh scan|annotate"; exit 1
        ;;
esac
