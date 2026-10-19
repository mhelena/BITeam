#!/usr/bin/env bash
# Merges QC'd study data with 1KGP3, excludes high-LD regions, LD prunes and applies per-ancestry HWE.
# Usage: bash 03_merge_and_prune.sh <study_prefix> <ref_prefix> <ref_samples> <out_prefix> <high_ld_file> <window_kb> <step> <r2> <maf> <threads> <ref_labels> [hwe_p]
set -euo pipefail

STUDY_PREFIX=$1
REF_PREFIX=$2           # PLINK1 or PLINK2 prefix
REF_PSAM=$3             # FID IID
OUT_PREFIX=$4
HIGH_LD_FILE=$5
WINDOW_KB=$6
STEP=$7
R2=$8
MAF=$9
THREADS=${10}
REF_LABELS=${11}        # IID, SuperPop
HWE_PVAL=${12:-1e-5}

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

awk '{print $1"\t"$4"\t"$5"\t"$6"\t"$2}' "${STUDY_PREFIX}.bim" \
    | sort -k1,1 -k2,2n > "${OUT_PREFIX}_study_snps.txt"

if [[ -f "${REF_PREFIX}.bim" ]]; then
    awk '{print $1"\t"$4"\t"$5"\t"$6"\t"$2}' "${REF_PREFIX}.bim" \
        | sort -k1,1 -k2,2n > "${OUT_PREFIX}_ref_snps.txt"
elif [[ -f "${REF_PREFIX}.pvar" ]]; then
    awk 'NR>1 && !/^#/{print $1"\t"$2"\t"$4"\t"$5"\t"$3}' "${REF_PREFIX}.pvar" \
        | sort -k1,1 -k2,2n > "${OUT_PREFIX}_ref_snps.txt"
fi

# Match on chr:pos:ref:alt in either allele order
python3 - <<PYEOF
import pandas as pd

study = pd.read_csv("${OUT_PREFIX}_study_snps.txt", sep='\t',
                    names=['chr','pos','ref','alt','snp_id'])
ref   = pd.read_csv("${OUT_PREFIX}_ref_snps.txt",   sep='\t',
                    names=['chr','pos','ref','alt','snp_id'])

def norm_chr(c):
    return str(c).replace('chr','')

study['chr'] = study['chr'].apply(norm_chr)
ref['chr']   = ref['chr'].apply(norm_chr)

m = pd.merge(study, ref, on=['chr','pos','ref','alt'], suffixes=('_s','_r'))

ref_flip = ref.copy()
ref_flip['ref'], ref_flip['alt'] = ref['alt'], ref['ref']
m2 = pd.merge(study, ref_flip, on=['chr','pos','ref','alt'], suffixes=('_s','_r'))

combined = pd.concat([m[['snp_id_s','snp_id_r']], m2[['snp_id_s','snp_id_r']]])
combined = combined.drop_duplicates(subset='snp_id_s')

print(f"Study SNPs: {len(study)}")
print(f"Reference SNPs: {len(ref)}")
print(f"Overlapping SNPs: {len(combined)}")

combined['snp_id_s'].to_csv("${OUT_PREFIX}_overlap_study_ids.txt", index=False, header=False)
combined['snp_id_r'].to_csv("${OUT_PREFIX}_overlap_ref_ids.txt",   index=False, header=False)
PYEOF

log "Overlapping SNPs: $(wc -l < ${OUT_PREFIX}_overlap_study_ids.txt)"

plink \
    --bfile "$STUDY_PREFIX" \
    --extract "${OUT_PREFIX}_overlap_study_ids.txt" \
    --make-bed \
    --keep-allele-order \
    --threads "$THREADS" \
    --out "${OUT_PREFIX}_study_overlap"

# MAF in 1KGP3 is applied here
if [[ -f "${REF_PREFIX}.bed" ]]; then
    plink \
        --bfile "$REF_PREFIX" \
        --extract "${OUT_PREFIX}_overlap_ref_ids.txt" \
        --keep "$REF_PSAM" \
        --maf "$MAF" \
        --make-bed \
        --keep-allele-order \
        --threads "$THREADS" \
        --out "${OUT_PREFIX}_ref_overlap"
elif [[ -f "${REF_PREFIX}.pgen" ]]; then
    plink2 \
        --pfile "$REF_PREFIX" \
        --extract "${OUT_PREFIX}_overlap_ref_ids.txt" \
        --keep "$REF_PSAM" \
        --autosome \
        --maf "$MAF" \
        --make-bed \
        --threads "$THREADS" \
        --out "${OUT_PREFIX}_ref_overlap"
fi

log "Harmonising study alleles to the reference"
python3 - <<PYEOF
import pandas as pd

study_bim = pd.read_csv("${OUT_PREFIX}_study_overlap.bim", sep='\t',
    names=['chr','snp','cm','pos','a1','a2'])
ref_bim = pd.read_csv("${OUT_PREFIX}_ref_overlap.bim", sep='\t',
    names=['chr','snp','cm','pos','a1','a2'])

ref_dict = {}
for _, r in ref_bim.iterrows():
    key = (str(r.chr).replace('chr',''), r.pos)
    ref_dict[key] = (r.a1, r.a2, r.snp)

flip_list = []
remove_list = []

for _, s in study_bim.iterrows():
    key = (str(s.chr).replace('chr',''), s.pos)
    if key not in ref_dict:
        continue
    ra1, ra2, rid = ref_dict[key]
    if s.a1 == ra1 and s.a2 == ra2:
        pass
    elif s.a1 == ra2 and s.a2 == ra1:
        flip_list.append(s.snp)
    else:
        remove_list.append(s.snp)

print(f"SNPs to flip: {len(flip_list)}")
print(f"SNPs to remove (allele mismatch): {len(remove_list)}")

with open("${OUT_PREFIX}_flip.txt", 'w') as f:
    f.write('\n'.join(flip_list) + '\n')
with open("${OUT_PREFIX}_mismatch.txt", 'w') as f:
    f.write('\n'.join(remove_list) + '\n')
PYEOF

plink \
    --bfile "${OUT_PREFIX}_study_overlap" \
    --flip "${OUT_PREFIX}_flip.txt" \
    --exclude "${OUT_PREFIX}_mismatch.txt" \
    --make-bed \
    --keep-allele-order \
    --threads "$THREADS" \
    --out "${OUT_PREFIX}_study_harmonized"

# If --bmerge fails on remaining strand problems, drop the .missnp SNPs from both sides and retry
plink \
    --bfile "${OUT_PREFIX}_study_harmonized" \
    --bmerge "${OUT_PREFIX}_ref_overlap" \
    --make-bed \
    --keep-allele-order \
    --threads "$THREADS" \
    --out "${OUT_PREFIX}_merged" || {
    log "Merge failed, retrying without .missnp SNPs"
    if [[ -f "${OUT_PREFIX}_merged-merge.missnp" ]]; then
        plink --bfile "${OUT_PREFIX}_study_harmonized" \
              --exclude "${OUT_PREFIX}_merged-merge.missnp" \
              --make-bed --keep-allele-order --threads "$THREADS" \
              --out "${OUT_PREFIX}_study_harmonized2"
        plink --bfile "${OUT_PREFIX}_ref_overlap" \
              --exclude "${OUT_PREFIX}_merged-merge.missnp" \
              --make-bed --keep-allele-order --threads "$THREADS" \
              --out "${OUT_PREFIX}_ref_overlap2"
        plink --bfile "${OUT_PREFIX}_study_harmonized2" \
              --bmerge "${OUT_PREFIX}_ref_overlap2" \
              --make-bed --keep-allele-order --threads "$THREADS" \
              --out "${OUT_PREFIX}_merged"
    fi
}

log "Merged SNPs: $(wc -l < ${OUT_PREFIX}_merged.bim)"

awk 'NF==4 && !/^#/' "$HIGH_LD_FILE" \
    | awk '{gsub("chr",""); print $1, $2, $3, $4}' \
    > "${OUT_PREFIX}_high_ld_plink.txt"

plink \
    --bfile "${OUT_PREFIX}_merged" \
    --exclude range "${OUT_PREFIX}_high_ld_plink.txt" \
    --make-bed \
    --keep-allele-order \
    --threads "$THREADS" \
    --out "${OUT_PREFIX}_nohighld"

log "SNPs after high-LD exclusion: $(wc -l < ${OUT_PREFIX}_nohighld.bim)"

# The prune list is computed on reference samples only, so the study cohort's
# composition does not influence which SNPs go into the PCA.
log "LD pruning on reference samples (r2=${R2}, window=${WINDOW_KB}kb, step=${STEP})"
plink \
    --bfile "${OUT_PREFIX}_nohighld" \
    --keep "$REF_PSAM" \
    --make-bed \
    --keep-allele-order \
    --threads "$THREADS" \
    --out "${OUT_PREFIX}_nohighld_refonly"

plink \
    --bfile "${OUT_PREFIX}_nohighld_refonly" \
    --indep-pairwise "${WINDOW_KB}" "${STEP}" "${R2}" \
    --threads "$THREADS" \
    --out "${OUT_PREFIX}_prune"

plink \
    --bfile "${OUT_PREFIX}_nohighld" \
    --extract "${OUT_PREFIX}_prune.prune.in" \
    --make-bed \
    --keep-allele-order \
    --threads "$THREADS" \
    --out "${OUT_PREFIX}_pruned"

log "SNPs after LD pruning: $(wc -l < ${OUT_PREFIX}_pruned.bim)"

# HWE is tested within AFR, EAS, EUR and SAS 1KGP3 samples separately, since a
# pooled test would flag SNPs that differ in frequency between groups (Wahlund
# effect). AMR is left out because admixture alone breaks HWE there. A SNP that
# fails in any group is removed.
python3 - "$REF_LABELS" "$OUT_PREFIX" <<'PYEOF'
import pandas as pd, sys

labels_file = sys.argv[1]
out_prefix  = sys.argv[2]

df = pd.read_csv(labels_file, sep='\t')

col_map = {}
for c in df.columns:
    cl = c.lower().replace(' ','_').replace('-','_').strip('#')
    if cl in ('iid', 'sample', 'sample_id', 'individual_id'):
        col_map[c] = 'IID'
    elif 'super' in cl and 'pop' in cl:
        col_map[c] = 'SuperPop'
df.rename(columns=col_map, inplace=True)

if 'IID' not in df.columns or 'SuperPop' not in df.columns:
    print(f"ERROR: need IID and SuperPop columns. Found: {list(df.columns)}")
    sys.exit(1)

for pop in ['AFR', 'EAS', 'EUR', 'SAS']:
    sub = df[df['SuperPop'] == pop][['IID']].copy()
    sub['FID'] = 0
    out = f"{out_prefix}_hwe_{pop}_samples.txt"
    sub[['FID','IID']].to_csv(out, sep='\t', index=False, header=False)
    print(f"{pop}: {len(sub)} samples")
PYEOF

HWE_FAIL_ALL="${OUT_PREFIX}_hwe_fail_all.txt"
> "$HWE_FAIL_ALL"

for POP in AFR EAS EUR SAS; do
    SAMP_FILE="${OUT_PREFIX}_hwe_${POP}_samples.txt"
    if [[ ! -f "$SAMP_FILE" ]] || [[ ! -s "$SAMP_FILE" ]]; then
        log "Skipping $POP, no samples"
        continue
    fi

    plink \
        --bfile "${OUT_PREFIX}_pruned" \
        --keep "$SAMP_FILE" \
        --hardy \
        --threads "$THREADS" \
        --out "${OUT_PREFIX}_hwe_${POP}" > /dev/null 2>&1

    # .hwe columns: CHR SNP TEST A1 A2 GENO O(HET) E(HET) P
    if [[ -f "${OUT_PREFIX}_hwe_${POP}.hwe" ]]; then
        awk -v pval="$HWE_PVAL" 'NR>1 && $9+0 < pval+0 {print $2}' \
            "${OUT_PREFIX}_hwe_${POP}.hwe" >> "$HWE_FAIL_ALL"
        N_FAIL=$(awk -v pval="$HWE_PVAL" 'NR>1 && $9+0 < pval+0' "${OUT_PREFIX}_hwe_${POP}.hwe" | wc -l)
        log "HWE failures in $POP ($(wc -l < "$SAMP_FILE") samples): $N_FAIL"
    fi
done

sort -u "$HWE_FAIL_ALL" > "${OUT_PREFIX}_hwe_fail_unique.txt"
N_HWE_FAIL=$(wc -l < "${OUT_PREFIX}_hwe_fail_unique.txt")
log "Unique SNPs failing HWE in any group: $N_HWE_FAIL"

if [[ "$N_HWE_FAIL" -gt 0 ]]; then
    plink \
        --bfile "${OUT_PREFIX}_pruned" \
        --exclude "${OUT_PREFIX}_hwe_fail_unique.txt" \
        --make-bed \
        --keep-allele-order \
        --threads "$THREADS" \
        --out "${OUT_PREFIX}_pruned_hwe"
else
    cp "${OUT_PREFIX}_pruned.bed" "${OUT_PREFIX}_pruned_hwe.bed"
    cp "${OUT_PREFIX}_pruned.bim" "${OUT_PREFIX}_pruned_hwe.bim"
    cp "${OUT_PREFIX}_pruned.fam" "${OUT_PREFIX}_pruned_hwe.fam"
fi

log "SNPs after per-ancestry HWE: $(wc -l < ${OUT_PREFIX}_pruned_hwe.bim) (${OUT_PREFIX}_pruned_hwe)"
