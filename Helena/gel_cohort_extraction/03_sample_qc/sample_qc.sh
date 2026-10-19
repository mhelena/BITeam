#!/usr/bin/env bash
#BSUB -q short
#BSUB -M 8000
#BSUB -n 2
#BSUB -R "rusage[mem=8000]"
# Removes Ti/Tv outlier samples and variants failing HWE (p < 1e-5) from
# final/${COHORT}_ALL. Submitted by submit.sh.
set -euo pipefail
: "${PIPELINE_CONFIG:?export PIPELINE_CONFIG=/absolute/path/to/config.sh first}"
source "$PIPELINE_CONFIG"

module load "$PLINK2_MODULE"

IN=${WORKDIR}/final/${COHORT}_ALL
QC=${WORKDIR}/qc
mkdir -p "$QC"

plink2 --bfile "$IN" \
       --sample-counts cols=fid,ts,tv \
       --threads 2 --memory 6000 \
       --out "${QC}/${COHORT}_titv"

# Palindromic SNPs were removed in 02_extract_plink, which pushes Ti/Tv up to ~4, so
# outliers are flagged at mean +/- 3 SD instead of the usual 1.9-2.3 window.
awk 'BEGIN{OFS="\t"; print "FID","IID","TS","TV","TiTv"}
     NR>1 {print $1, $2, $3, $4, ($4>0 ? $3/$4 : "NA")}' \
    "${QC}/${COHORT}_titv.scount" > "${QC}/${COHORT}_titv.tsv"

read -r N MEAN SD < <(awk 'NR>1 && $5!="NA" {n++; s+=$5; ss+=$5*$5}
    END {m=s/n; printf "%d %.4f %.4f\n", n, m, sqrt(ss/n-m*m)}' "${QC}/${COHORT}_titv.tsv")
LO=$(awk -v m="$MEAN" -v s="$SD" 'BEGIN{printf "%.4f", m - 3*s}')
HI=$(awk -v m="$MEAN" -v s="$SD" 'BEGIN{printf "%.4f", m + 3*s}')

awk -v lo="$LO" -v hi="$HI" 'BEGIN{OFS="\t"}
    NR>1 && $5!="NA" && ($5<lo || $5>hi) {print $1, $2}' \
    "${QC}/${COHORT}_titv.tsv" > "${QC}/${COHORT}_titv_outliers.txt"
N_FAIL=$(wc -l < "${QC}/${COHORT}_titv_outliers.txt" | tr -d ' ')

cat > "${QC}/${COHORT}_titv_summary.txt" <<EOF
Per-sample Ti/Tv, ${COHORT}_ALL
  n_samples : ${N}
  mean Ti/Tv: ${MEAN}
  SD Ti/Tv  : ${SD}
  threshold : ${LO} to ${HI} (mean +/- 3 SD)
  outliers  : ${N_FAIL}
EOF
cat "${QC}/${COHORT}_titv_summary.txt"

REMOVE_OPT=""
if [[ "$N_FAIL" -gt 0 ]]; then
    REMOVE_OPT="--remove ${QC}/${COHORT}_titv_outliers.txt"
fi

# shellcheck disable=SC2086
plink2 --bfile "$IN" \
       $REMOVE_OPT \
       --hwe 1e-5 \
       --hardy \
       --make-bed \
       --threads 2 --memory 6000 \
       --out "${QC}/${COHORT}_qc"

N_IN_VAR=$(wc -l < "${IN}.bim" | tr -d ' ')
N_IN_SAM=$(wc -l < "${IN}.fam" | tr -d ' ')
N_OUT_VAR=$(wc -l < "${QC}/${COHORT}_qc.bim" | tr -d ' ')
N_OUT_SAM=$(wc -l < "${QC}/${COHORT}_qc.fam" | tr -d ' ')

cat > "${QC}/${COHORT}_qc_summary.txt" <<EOF
${COHORT} sample QC
  variants: ${N_IN_VAR} -> ${N_OUT_VAR} (HWE p < 1e-5 removed $((N_IN_VAR - N_OUT_VAR)))
  samples : ${N_IN_SAM} -> ${N_OUT_SAM} (Ti/Tv removed $((N_IN_SAM - N_OUT_SAM)))
EOF
cat "${QC}/${COHORT}_qc_summary.txt"
