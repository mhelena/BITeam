#!/usr/bin/env bash
# Settings for the GEL RE burden-VCF build. burden_submit.sh reads the config.sh
# next to it, or the file in $PIPELINE_CONFIG if that is set.

export LSF_PROJECT=re_gecip_neurology

export PROJECT_DIR="/re_gecip/<your_gecip>/<project>"
[[ "$PROJECT_DIR" != *"<"* ]] || { echo "config.sh: set PROJECT_DIR first" >&2; return 1; }

export OUTDIR="${PROJECT_DIR}/burden_vcf"

# Cohorts to convert and merge, in the sample order wanted in the merged VCF.
export COHORTS="dee epi_controls healthy_controls"

# Directory holding <cohort>_chr1..22.{bed,bim,fam} for one cohort. The default
# matches the layout written by gel_cohort_extraction.
perchr_dir() { echo "${PROJECT_DIR}/$1/scratch/perchr"; }

# Merged VCF is ${OUTDIR}/${MERGE_OUT_TAG}_for_burden.vcf.gz
export MERGE_OUT_TAG="controls_3cohort"

# Input of the acgt0 step: the merged VCF after variant QC.
export ACGT0_IN="${OUTDIR}/${MERGE_OUT_TAG}_for_burden_postQC.vcf.gz"

# Optional: copy the acgt0 output here and check the copy (md5 and record count).
export OUTPUT_COPY_DIR=""

# The acgt0 job needs a queue allowing at least 12 h; check with `bqueues -u $USER`.
export ACGT0_QUEUE=medium
