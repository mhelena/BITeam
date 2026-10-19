#!/usr/bin/env bash
# Settings read by every script. Edit the paths, then point the scripts at it with
#   export PIPELINE_CONFIG=/absolute/path/to/config.sh

export LSF_PROJECT=re_gecip_neurology

export PROJECT_DIR="/re_gecip/<your_gecip>/<project>"
[[ "$PROJECT_DIR" != *"<"* ]] || { echo "config.sh: set PROJECT_DIR first" >&2; return 1; }

# Folder holding these scripts on the RE (defaults to the folder of this file).
export SCRIPT_DIR=${SCRIPT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}

# PanelApp panel (402 = Early onset or syndromic epilepsy) and the confidence
# level kept: Green, Amber (Green + Amber) or anything else for all levels.
export PANEL_ID=402
export PANEL_VERSION=8.167
export PANEL_CONFIDENCE=Green

# How the scan picks panel variants:
#   bed     restrict to PANEL_BED with bcftools -R (any gene annotated there)
#   symbol  scan whole chunks and keep rows whose VEP SYMBOL is in PANEL_GENES
export PANEL_MODE=bed

export WORK_DIR=${PROJECT_DIR}/panel_${PANEL_ID}_v${PANEL_VERSION}
export PANEL_GENES=${WORK_DIR}/panel_genes.txt
export PANEL_BED=${WORK_DIR}/panel_cds_hg38.bed
export SITES_DIR=${WORK_DIR}/sites
export ANNOT_DIR=${WORK_DIR}/annot
export FINAL_DIR=${WORK_DIR}/final
export LOG_DIR=${WORK_DIR}/logs

# Candidate controls, one platekey per line. Carriers are removed from this list.
export POOL_FILE=${PROJECT_DIR}/controls/pool_candidates.txt
export FINAL_OUT=${FINAL_DIR}/controls.txt

# Optional, for panel_3_final.R only: the lists that were unioned into
# POOL_FILE (space-separated), so survival is reported per list. Platekeys or
# participant ids; ids are mapped to platekeys through LabKey.
export POOL_SUBSETS="${POOL_FILE}"
export LABKEY_BASE_URL=https://labkey.prod.aws.gel.ac/labkey
export LABKEY_FOLDER=/main-programme/main-programme_v19_2024-10-31

export AGGV2=/gel_data_resources/main_programme/aggregation/aggregate_gVCF_strelka/aggV2
export FA_DIR=${AGGV2}/functional_annotation/VEP_109
export GENO_DIR=${AGGV2}/genomic_data
export GENCODE_GTF=/public_data_resources/gencode/v38/gencode.v38.primary_assembly.annotation.gtf.gz

export BCFTOOLS_MODULE=bcftools/1.16
# bcftools/1.16 on the RE does not ship bgzip or tabix
export HTSLIB_MODULE=htslib/1.21
export BEDTOOLS_MODULE=bedtools/2.30.0
