#!/usr/bin/env bash
# Settings for the SGE runner. qc_submit.sh reads the config.sh next to it, or the
# file in $PIPELINE_CONFIG if that is set.

# allanc: all samples, MEDIAN_DP/MEDIAN_GQ filter, no HWE.
# eur:    samples in EUR_KEEP only, HWE p < 1e-5, no MEDIAN_DP/MEDIAN_GQ filter.
# A MODE already set in the environment wins: MODE=eur ./qc_submit.sh
export MODE=${MODE:-allanc}

export WORKDIR="/SAN/<project>/<user>/<study>"
[[ "$WORKDIR" != *"<"* ]] || { echo "config.sh: set WORKDIR first" >&2; return 1; }

# Joint-called, PASS-only, GRCh38 multi-sample VCF (bgzipped and indexed).
export RAW_VCF="/SAN/<project>/<path>/cohort.sorted.PASS_only.vcf.gz"
export OUTDIR="${WORKDIR}/vcf_${MODE}"

# eur mode only. PLINK keep file (FID IID); the IID column is used.
export EUR_KEEP="${WORKDIR}/eur_keep.txt"

# Optional repeat/low-complexity mask, tabix-indexed BED, e.g. GIAB
# GRCh38_AllTandemRepeatsandHomopolymers_slop5.bed.gz. Empty skips the mask.
export LCR_BED=""

# Sourced before submitting (qsub -V passes it on to the jobs); sets
# LD_LIBRARY_PATH for the shared python3 build on UCL CS.
export PYTHON_ENV_SOURCE=/share/apps/source_files/python/python-3.12.11.source

# There is no module system on the UCL CS cluster, so tools need full paths.
# bcftools 1.9 is enough; HWE is computed by helpers/hwe_exact.py for that reason.
export BCFTOOLS="${BCFTOOLS:-/cluster/project8/<group>/<path>/bcftools-1.9/bcftools}"
export BGZIP="${BGZIP:-/cluster/project8/SMS_Epilepsy_Genetics/Software/samtools-1.2/htslib-1.2.1/bgzip}"
export TABIX="${TABIX:-/cluster/project8/SMS_Epilepsy_Genetics/Software/samtools-1.2/htslib-1.2.1/tabix}"
export PYTHON3="${PYTHON3:-python3}"
