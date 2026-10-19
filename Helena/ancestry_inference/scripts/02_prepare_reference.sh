#!/usr/bin/env bash
# Downloads the 1KGP3 reference panel and writes the sample list and super-population labels.
# Usage: bash 02_prepare_reference.sh <out_dir> <hg38|hg19> <threads>
set -euo pipefail

OUT_DIR=$1
GENOME_BUILD=$2
THREADS=$3

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }
mkdir -p "$OUT_DIR"

if [[ "$GENOME_BUILD" == "hg38" ]]; then
    # PLINK2 resource-page copies of 1KGP3 on GRCh38
    PLINK2_URL="https://www.dropbox.com/s/j72j6uciq5zuzii/all_hg38.pgen.zst?dl=1"
    PVAR_URL="https://www.dropbox.com/scl/fi/id642dpdd858uy41og8qi/all_hg38_rs_noannot.pvar.zst?rlkey=sskyiyam1bsqweujjmxqv1h55&dl=1"
    PSAM_URL="https://www.dropbox.com/scl/fi/u5udzzaibgyvxzfnjcvjc/hg38_corrected.psam?rlkey=oecjnk4vmbhc8b1p202l0ih4x&dl=1"
    log "Downloading 1KGP3 hg38 PLINK2 files"
    wget -q -O "${OUT_DIR}/all_hg38.pgen.zst"  "$PLINK2_URL" || \
        { log "ERROR: download failed, fetch 1KGP3 manually and set ref_1kg_plink"; exit 1; }
    wget -q -O "${OUT_DIR}/all_hg38.pvar.zst"  "$PVAR_URL"
    wget -q -O "${OUT_DIR}/hg38_corrected.psam" "$PSAM_URL"
    plink2 --zst-decompress "${OUT_DIR}/all_hg38.pgen.zst" > "${OUT_DIR}/all_hg38.pgen"
else
    log "Downloading 1KGP3 hg19 VCFs"
    for CHR in $(seq 1 22); do
        URL="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
        wget -q -O "${OUT_DIR}/1kg_chr${CHR}.vcf.gz" "$URL" || true
    done
    wget -q -O "${OUT_DIR}/integrated_call_samples_v3.20130502.ALL.panel" \
        "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel"
fi

log "Writing sample list and super-population labels"
python3 - "$OUT_DIR" <<'PYEOF'
import pandas as pd
import sys, os

out_dir = sys.argv[1]

panel_files = [
    f"{out_dir}/hg38_corrected.psam",
    f"{out_dir}/integrated_call_samples_v3.20130502.ALL.panel",
]
df = None
for f in panel_files:
    if os.path.exists(f):
        df = pd.read_csv(f, sep='\t')
        print(f"Loaded: {f}")
        break

if df is None:
    print("No panel file found")
    sys.exit(0)

rename = {}
for c in df.columns:
    cl = c.strip('#').lower()
    if cl in ('iid', 'sample', 'sample_name', 'sample_id', 'individual_id'):
        rename[c] = 'IID'
    elif cl in ('super_pop', 'superpop', 'super population code'):
        rename[c] = 'SuperPop'
    elif cl in ('pop', 'population', 'population_code'):
        rename[c] = 'SubPop'
df.rename(columns=rename, inplace=True)

if 'IID' not in df.columns:
    df.columns = [c.strip('#') for c in df.columns]
    if '#IID' in df.columns: df.rename(columns={'#IID': 'IID'}, inplace=True)

df['FID'] = df['IID']
df[['FID','IID']].to_csv(f"{out_dir}/1kg_unrelated_samples.txt", sep='\t', index=False, header=False)
print(f"Written {len(df)} samples")
if 'SuperPop' in df.columns:
    df[['IID','SuperPop']].to_csv(f"{out_dir}/1kg_super_pop_labels.tsv", sep='\t', index=False)
PYEOF

log "Reference data in $OUT_DIR"
