# Settings for the gene carrier screen. All scripts read this file: export
# PIPELINE_CONFIG=/absolute/path/to/config.sh before running or submitting anything.

export LSF_PROJECT="re_gecip_neurology"
export LSF_QUEUE="short"

# Gene of interest, GRCh38. REGION should cover the whole gene with a small margin.
export GENE="SCN1A"
export REGION="chr2:165990000-166150000"
# Optional Ensembl transcript ID to pin the annotation to (e.g. the MANE Select
# transcript). Empty means use the VEP CANONICAL flag.
export TRANSCRIPT=""

# Participants matching the phenotype (HPO regex or GEL disease_sub_group, several
# subgroups separated by ";") are removed from the candidate pool.
export PHENO_HPO_REGEX="epilep|seizure|convuls"
export PHENO_DISEASE_SUBGROUP="Inherited Epilepsy Syndromes"

# Disease categories unrelated to the phenotype, one per line, matched against GEL
# specific_disease, disease_sub_group, disease_group and HPO term. Only carriers
# matching this list are kept. For epilepsy: Supplementary Table 3 of PMID 37006128.
export WHITELIST="/re_gecip/YOUR_GECIP/YOUR_USER/gene_carrier_screen/non_epi_diseases.txt"

# Optional file of platekeys or participant IDs to drop, one per line (e.g. carriers
# removed after manual review). Leave empty to keep everyone.
export EXCLUDE_IDS=""

# Age filter for the HPO table and review sheet: at least MIN_AGE at the data cut.
export DATA_CUT_YEAR=2024
export MIN_AGE=4

export LABKEY_BASE="https://labkey.prod.aws.gel.ac/labkey"
export LABKEY_FOLDER="/main-programme/main-programme_v19_2024-10-31"

AGGV2="/gel_data_resources/main_programme/aggregation/aggregate_gVCF_strelka/aggV2"
export FA_DIR="${AGGV2}/functional_annotation/VEP_109"
export GENO_DIR="${AGGV2}/genomic_data"

export OUT_DIR="/re_gecip/YOUR_GECIP/YOUR_USER/gene_carrier_screen/${GENE}"
export WORK="${OUT_DIR}/work"
export LOG_DIR="${OUT_DIR}/logs"

# Files passed between steps; normally no need to change these.
export POOL_FILE="${OUT_DIR}/pool.txt"
export CARRIERS="${OUT_DIR}/carriers.txt"
export CARRIER_VARIANTS="${OUT_DIR}/carrier_variants.tsv"
export VARIANTS="${WORK}/variants.tsv"
export VAR_ANNOT="${OUT_DIR}/variants_annotated.tsv"
export CONTROLS_LIST="${OUT_DIR}/carriers_whitelist_eur08.txt"
export HPO_WIDE="${OUT_DIR}/carriers_whitelist_eur08_hpo_wide.tsv"
# Sample list for the Airlock table. Point this at the post-review list if it differs.
export FINAL_SAMPLES="${CONTROLS_LIST}"
export AIRLOCK_DIR="${OUT_DIR}/for_airlock"

REGION_CHR="${REGION%%:*}"
_range="${REGION#*:}"
export REGION_CHR REGION_START="${_range%-*}" REGION_END="${_range#*-}"
unset _range
