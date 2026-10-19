#!/usr/bin/env bash
# PCA on unrelated 1KGP3 samples, then projects study samples onto the reference PCs.
# Usage: bash 04_pca_gcta.sh <merged_prefix> <ref_samples> <study_samples> <out_prefix> <n_pcs> <king_cutoff> <threads>
# gcta64 here is bin/gcta64, a plink2 wrapper, so bin/ must be on PATH.
set -euo pipefail

MERGED_PREFIX=$1
REF_SAMPLES=$2
STUDY_SAMPLES=$3
OUT_PREFIX=$4
N_PCS=$5
KING_CUTOFF=$6
THREADS=$7

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

# The 1KGP3 release is meant to be unrelated already; KING-robust is a check.
log "KING-robust kinship on reference samples"
plink2 \
    --bfile "$MERGED_PREFIX" \
    --keep "$REF_SAMPLES" \
    --make-bed \
    --threads "$THREADS" \
    --out "${OUT_PREFIX}_ref_only"

plink2 \
    --bfile "${OUT_PREFIX}_ref_only" \
    --make-king-table \
    --king-table-filter "$KING_CUTOFF" \
    --threads "$THREADS" \
    --out "${OUT_PREFIX}_ref_king"

if [[ -f "${OUT_PREFIX}_ref_king.kin0" ]]; then
    python3 - <<PYEOF
import pandas as pd, os
from collections import Counter

kin_file = "${OUT_PREFIX}_ref_king.kin0"
ref_file = "$REF_SAMPLES"

if os.path.getsize(kin_file) == 0:
    print("No related pairs")
    import shutil; shutil.copy(ref_file, "${OUT_PREFIX}_ref_unrelated.txt")
else:
    kin = pd.read_csv(kin_file, sep='\t')
    # Greedy: from each related pair drop the sample with more relatives
    pairs = list(zip(kin['#FID1'].astype(str)+'_'+kin['IID1'].astype(str),
                     kin['FID2'].astype(str)+'_'+kin['IID2'].astype(str)))
    counts = Counter(s for pair in pairs for s in pair)
    remove = set()
    for s1, s2 in pairs:
        if s1 in remove or s2 in remove:
            continue
        if counts[s1] >= counts[s2]:
            remove.add(s1)
        else:
            remove.add(s2)
    print(f"Removing {len(remove)} related reference samples")
    ref = pd.read_csv(ref_file, sep='\t', names=['FID','IID'])
    ref['key'] = ref['FID'].astype(str)+'_'+ref['IID'].astype(str)
    ref_unrel = ref[~ref['key'].isin(remove)][['FID','IID']]
    ref_unrel.to_csv("${OUT_PREFIX}_ref_unrelated.txt", sep='\t', index=False, header=False)
    print(f"Unrelated reference samples: {len(ref_unrel)}")
PYEOF
else
    cp "$REF_SAMPLES" "${OUT_PREFIX}_ref_unrelated.txt"
fi

log "PCA on ${N_PCS} PCs"
plink \
    --bfile "$MERGED_PREFIX" \
    --keep "${OUT_PREFIX}_ref_unrelated.txt" \
    --make-bed \
    --keep-allele-order \
    --threads "$THREADS" \
    --out "${OUT_PREFIX}_ref_for_pca"

gcta64 \
    --bfile "${OUT_PREFIX}_ref_for_pca" \
    --pca "$N_PCS" \
    --thread-num "$THREADS" \
    --out "${OUT_PREFIX}_ref_pca"

# Projecting rather than running a joint PCA keeps the PC axes defined by the
# reference alone, which is what the classifier is trained on.
log "Projecting study samples"
plink \
    --bfile "$MERGED_PREFIX" \
    --keep "$STUDY_SAMPLES" \
    --make-bed \
    --keep-allele-order \
    --threads "$THREADS" \
    --out "${OUT_PREFIX}_study_for_proj"

gcta64 \
    --bfile "${OUT_PREFIX}_study_for_proj" \
    --project-loading "${OUT_PREFIX}_ref_pca" "$N_PCS" \
    --thread-num "$THREADS" \
    --out "${OUT_PREFIX}_study_projected"

python3 - <<PYEOF
import pandas as pd

ref_pcs = pd.read_csv("${OUT_PREFIX}_ref_pca.eigenvec", sep=r'\s+', header=None)
n_pcs = ref_pcs.shape[1] - 2
ref_pcs.columns = ['FID','IID'] + [f'PC{i}' for i in range(1, n_pcs+1)]
ref_pcs['source'] = 'reference_1KGP3'

study_pcs = pd.read_csv("${OUT_PREFIX}_study_projected.proj.eigenvec", sep=r'\s+', header=None)
study_pcs.columns = ['FID','IID'] + [f'PC{i}' for i in range(1, n_pcs+1)]
study_pcs['source'] = 'study'

all_pcs = pd.concat([ref_pcs, study_pcs], ignore_index=True)
all_pcs.to_csv("${OUT_PREFIX}_all_pcs.tsv", sep='\t', index=False)
print(f"PCs: {len(ref_pcs)} reference + {len(study_pcs)} study")
PYEOF

log "PCs written to ${OUT_PREFIX}_all_pcs.tsv"
