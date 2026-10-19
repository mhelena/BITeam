#!/usr/bin/env bash
# Settings read by every script. Edit the paths, then point the scripts at it with
#   export PIPELINE_CONFIG=/absolute/path/to/config.sh

export LSF_PROJECT=re_gecip_neurology

export PROJECT_DIR="/re_gecip/<your_gecip>/<project>"
[[ "$PROJECT_DIR" != *"<"* ]] || { echo "config.sh: set PROJECT_DIR first" >&2; return 1; }

# A COHORT already set in the environment wins, so several cohorts can share
# this file: COHORT=dee ./02_extract_plink/submit.sh
export COHORT=${COHORT:-my_cohort}
export WORKDIR=${PROJECT_DIR}/${COHORT}

# Per-chunk intermediates, roughly 10-20 GB for a few hundred samples. Can point
# to /re_scratch instead. Safe to delete once final/${COHORT}_ALL.bed exists.
export SCRATCH=${WORKDIR}/scratch

export AGGV2_DIR=/gel_data_resources/main_programme/aggregation/aggregate_gVCF_strelka/aggV2/genomic_data

# One aggV2 chunk name per line (first column), e.g. chr1_<start>_<end>, such
# that ${AGGV2_DIR}/gel_mainProgramme_aggV2_<chunk>.vcf.gz exists.
export CHUNKS_TXT=${PROJECT_DIR}/reference/chunks.txt

# Tab-separated with header: platekey, participant_id, pred_european_ancestries, sex
export ANCESTRY_SEX_TSV=${PROJECT_DIR}/reference/aggV2_ancestry_sex.tsv

export BCFTOOLS_MODULE=bcftools/1.16
export PLINK2_MODULE=plink/2.0
