#!/usr/bin/env bash
# Settings for qc_local_run.sh. It reads the config.sh next to it, or the file in
# $PIPELINE_CONFIG if that is set.

# Joint-called, PASS-only, GRCh38 multi-sample VCF (bgzipped and indexed).
IN="/path/to/cohort.sorted.PASS_only.vcf.gz"
[[ "$IN" != /path/to/* ]] || { echo "config.sh: set IN first" >&2; return 1; }

OUTDIR="/path/to/results/qc_allanc"
OUT_PREFIX="cohort_allanc"

BCFTOOLS=${BCFTOOLS:-bcftools}
BGZIP=${BGZIP:-bgzip}
TABIX=${TABIX:-tabix}
PYTHON3=${PYTHON3:-python3}
