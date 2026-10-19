#!/bin/bash
# Runs the pipeline with a config profile (local or cluster). Run from this directory.
# Usage: bash run_pipeline.sh [local|cluster] [extra snakemake args], or on SGE: qsub run_pipeline.sh cluster
#$ -N ancestry_pipeline
#$ -cwd
#$ -pe smp 8
#$ -l h_vmem=4G
#$ -l h_rt=12:00:00
#$ -o logs/pipeline.log
#$ -e logs/pipeline.err

set -eo pipefail

PROFILE=${1:-local}
[[ $# -gt 0 ]] && shift

case "$PROFILE" in
    local)   CONFIG=config/config.yaml ;;
    cluster) CONFIG=config/cluster.yaml ;;
    *) echo "Unknown profile '$PROFILE' (use local or cluster)" >&2; exit 1 ;;
esac

# SGE jobs start in a non-interactive shell, so conda has to be sourced explicitly.
# e.g. CONDA_SH=/share/apps/miniconda-4.12.0/etc/profile.d/conda.sh CONDA_ENV=ancestry_inference
if [[ -n "${CONDA_SH:-}" ]]; then
    source "$CONDA_SH"
fi
if [[ -n "${CONDA_ENV:-}" ]]; then
    conda activate "$CONDA_ENV"
fi

mkdir -p logs

# bin/ must be on PATH for the gcta64 wrapper; plink/plink2/bcftools come from conda.
PATH="${CONDA_PREFIX:+$CONDA_PREFIX/bin:}$PWD/bin:$PATH" \
    snakemake -j "${NSLOTS:-8}" --rerun-incomplete --configfile "$CONFIG" "$@"
