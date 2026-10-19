#!/usr/bin/env bash
# Settings read by every script. Edit the paths, then point the scripts at it with
#   export PIPELINE_CONFIG=/absolute/path/to/config.sh

export LSF_PROJECT=re_gecip_neurology
export PLINK2_MODULE=plink/2.0
export PLINK19_MODULE=plink/1.9
export R_MODULE=R/4.3.3

export PROJECT_DIR="/re_gecip/<your_gecip>/<project>"
[[ "$PROJECT_DIR" != *"<"* ]] || { echo "config.sh: set PROJECT_DIR first" >&2; return 1; }

export COHORT=${COHORT:-my_cohort}
export WORKDIR=${PROJECT_DIR}/${COHORT}

# Cases after sample QC (output of gel_cohort_extraction), and the three control
# sets with cohort-level QC done. Joint QC happens in 02_target_qc.sh.
export CASE_PLINK=${WORKDIR}/qc/${COHORT}_qc
export HC_PLINK=${PROJECT_DIR}/controls/hc_new_qc
export DEE_PLINK=${PROJECT_DIR}/controls/dee_ALL_hwe
export EPI_PLINK=${PROJECT_DIR}/controls/epilepsy_controls

# Label for the cases in tables and plots
export CASE_LABEL=Cases

# Tab-separated with header: platekey, participant_id, pred_european_ancestries, sex
# (sex 1 male, 2 female). aggV2 .fam files carry sex 0, so sex comes from here.
export ANCESTRY_SEX_TSV=${PROJECT_DIR}/reference/aggV2_ancestry_sex.tsv

export KG_PLINK=${PROJECT_DIR}/reference/Merge2_forAncestry
export KG_PANEL=/public_data_resources/1000-genomes/integrated_call_samples_v3.20130502.ALL.panel

# Control cohort PRSice compares the cases against to pick the best P-value
# threshold. The other two are still scored and plotted. One of: epi, hc, dee
export PRSICE_CONTROL=epi

export PRSICE=${PROJECT_DIR}/PRSice/PRSice_linux
export BASE_HG38=${PROJECT_DIR}/base_hg38

# Scores to run, and one GWAS_<name> base file per score.
export PRS_LIST="epilepsy intelligence mvage"
export GWAS_epilepsy=${BASE_HG38}/Epilepsy.PRSice.hg38.tsv.gz
export GWAS_intelligence=${BASE_HG38}/Intelligence.PRSice.hg38.txt
export GWAS_longevity=${BASE_HG38}/Longevity.PRSice.hg38.tsv.gz
export GWAS_bmi=${BASE_HG38}/BMI.PRSice.hg38.txt
# The raw mvAge file has rsIDs and matches nothing in aggV2 (chr:pos:ref:alt IDs).
# If Num_SNP in the .summary is tiny, try mvAge.PRSice2.EUR.hg38.chrposID_alt.tsv
# (A1/A2 flipped).
export GWAS_mvage=${BASE_HG38}/mvAge.PRSice2.EUR.hg38.chrposID.tsv

# modules/apoe: which entry of PRS_LIST is the longevity score
export LONGEVITY_SCORE=longevity

# modules/within_case
export WITHIN_DIR=${WORKDIR}/within_case
# Binary outcome measured in the cases. OUTCOME_TSV is tab-separated with header,
# IID plus a column named like OUTCOME (any case) holding yes/no.
export OUTCOME=outcome
export OUTCOME_TSV=${WORKDIR}/phenotypes/outcome.tsv
# Optional continuous outcome: tab-separated, IID and onset_months. Leave empty to skip.
export ONSET_TSV=
# Tab-separated with header IID and Sex (male/female or 1/2)
export WITHIN_SEX_TSV=${WORKDIR}/phenotypes/sex.tsv
export WITHIN_PRS_LIST="bmi epilepsy intelligence longevity mvage"
# Scores that enter the joint within-case models (must be in WITHIN_PRS_LIST)
export JOINT_PRS_LIST="epilepsy intelligence mvage bmi"
# Optional: tab-separated, IID plus extra covariates (age, centre, batch...)
export EXTRA_COV_TSV=${WITHIN_DIR}/extra_cov.txt
export N_REP=100
export N_PERM=1000
