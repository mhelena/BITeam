#!/usr/bin/env bash
# Join the annotation rows with the pool allele counts on CHROM:POS:REF:ALT into
# _merged/panel_annot_with_pool_ac.tsv. Submitted by panel_annot_submit.sh.
#BSUB -q short
#BSUB -M 4000
#BSUB -n 1
#BSUB -R "rusage[mem=4000]"
set -euo pipefail
source "${PIPELINE_CONFIG:?set PIPELINE_CONFIG to the absolute path of config.sh}"

WORK="${ANNOT_DIR}/_merged"
annot="${WORK}/panel_annot_all.tsv"
ac="${WORK}/panel_pool_ac.tsv"
out="${WORK}/panel_annot_with_pool_ac.tsv"

[[ -s "$annot" ]] || { echo "missing $annot"; exit 1; }
[[ -s "$ac"    ]] || { echo "missing $ac";    exit 1; }

awk -F'\t' 'BEGIN { OFS="\t" }
    FNR == NR {
        key = $1 ":" $2 ":" $3 ":" $4
        ac[key]    = $5
        an[key]    = $6
        af[key]    = $7
        het[key]   = $8
        hom[key]   = $9
        next
    }
    {
        key = $1 ":" $2 ":" $3 ":" $4
        a = (key in ac) ? ac[key]  : "."
        b = (key in an) ? an[key]  : "."
        c = (key in af) ? af[key]  : "."
        d = (key in het) ? het[key]: "."
        e = (key in hom) ? hom[key]: "."
        print $0, a, b, c, d, e
    }
' "$ac" "$annot" > "$out"

hdr="CHROM\tPOS\tREF\tALT\tSYMBOL\tConsequence\tIMPACT\tgnomADg_AF\tPOOL_AC\tPOOL_AN\tPOOL_AF\tPOOL_NHET\tPOOL_NHOM"
{ printf "%b\n" "$hdr"; cat "$out"; } > "${out}.tmp" && mv "${out}.tmp" "$out"

n=$(($(wc -l < "$out") - 1))
echo "[$(date)] merged annotation: ${n} rows"
echo "         written to: $out"
