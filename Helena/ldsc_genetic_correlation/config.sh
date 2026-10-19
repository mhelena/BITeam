# Settings for the LDSC genetic correlation pipeline.
# Copy this file, edit the paths, and point PIPELINE_CONFIG at your copy.

LDSC_DIR=/path/to/ldsc                 # clone of github.com/bulik/ldsc
PYTHON=python3                         # needs to be the python env LDSC runs in

# European LD scores and regression weights (HapMap3, hg19)
REF_LD=$LDSC_DIR/LDscore/LDscore
W_LD=$LDSC_DIR/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC

# Optional: restrict munged SNPs to the HapMap3 list. Leave empty to keep all.
MERGE_ALLELES=

GWAS_DIR=/path/to/gwas_sumstats        # 'file' column in the trait table is relative to this
TRAITS=$(dirname "${BASH_SOURCE[0]}")/traits/traits_ds.tsv
OUTDIR=/path/to/ldsc_out

# Traits whose pairs go in the forest plot (comma separated names). Empty = all pairs.
PLOT_FOCUS=
# Colour scale limit for the heatmap
PLOT_RG_LIMIT=1
