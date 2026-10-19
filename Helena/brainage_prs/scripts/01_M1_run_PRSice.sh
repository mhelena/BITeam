#!/bin/bash
# Runs PRSice-2 for 8 GWAS bases x 4 binary contrasts (32 runs), writing .all_score and .summary per run.
# Usage: bash scripts/01_M1_run_PRSice.sh config.sh

CONFIG=${1:-$PIPELINE_CONFIG}
[ -f "$CONFIG" ] || { echo "usage: $0 config.sh (or set PIPELINE_CONFIG)" >&2; exit 1; }
source "$CONFIG"

for v in PRSICE TARGET VALID_SNPS COVARIATES PHENO_DIR PRS_DIR; do
  [ -n "${!v}" ] || { echo "missing $v in $CONFIG" >&2; exit 1; }
done
mkdir -p "$PRS_DIR"

N=0; PASS=0; FAIL=0
for gwas in epilepsy intelligence longevity AD brainage_gm brainage_wm brainage_gwm frailty; do
  base_var="GWAS_${gwas}"
  BASE="${!base_var}"
  [ -n "$BASE" ] || { echo "missing $base_var in $CONFIG" >&2; exit 1; }
  for pheno in syndrome DRE convulsions currentVPA; do
    N=$((N+1))
    OUT="$PRS_DIR/PRS_${gwas}_${pheno}"
    echo "[$N/32] $gwas x $pheno"
    "$PRSICE" \
      --A1 A1 --A2 A2 --bp BP --chr CHR \
      --snp SNP --stat BETA --pvalue P \
      --all-score --fastscore \
      --bar-levels 0.0001,0.001,0.01,0.05,0.1,0.5,1 \
      --base "$BASE" \
      --target "$TARGET" \
      --extract "$VALID_SNPS" \
      --cov-file "$COVARIATES" \
      --pheno-file "$PHENO_DIR/pheno_${pheno}.txt" \
      --pheno-col Pheno \
      --binary-target T \
      --out "$OUT" \
      --print-snp
    [ $? -eq 0 ] && PASS=$((PASS+1)) || FAIL=$((FAIL+1))
  done
done
echo "Done: $PASS ok / $FAIL failed (total $N)"
