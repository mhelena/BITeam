#!/usr/bin/env bash
# Variant QC of the study VCF (GEL AggV2 SNP QC; high-LD exclusion, LD pruning and HWE are in step 3).
# Usage: bash 01_qc_vcf.sh <in.vcf.gz> <out_prefix> <hg38|hg19> <maf> <max_missing> <min_gq> <min_dp> <threads>
set -euo pipefail

INPUT_VCF=$1
OUT_PREFIX=$2
GENOME_BUILD=$3
MAF=$4
MAX_MISS=$5
MIN_GQ=$6
MIN_DP=$7
THREADS=$8

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

if [[ "$GENOME_BUILD" == "hg38" ]]; then
    REGIONS=$(seq 1 22 | sed 's/^/chr/' | tr '\n' ',' | sed 's/,$//')
else
    REGIONS=$(seq 1 22 | tr '\n' ',' | sed 's/,$//')
fi

log "Keeping autosomal biallelic SNPs"
bcftools view \
    --regions "$REGIONS" \
    --type snps \
    --min-alleles 2 \
    --max-alleles 2 \
    "$INPUT_VCF" \
    -O z -o "${OUT_PREFIX}_s1_biallelic.vcf.gz"

bcftools index --tbi "${OUT_PREFIX}_s1_biallelic.vcf.gz"
log "Variants: $(bcftools view -H ${OUT_PREFIX}_s1_biallelic.vcf.gz | wc -l)"

# GEL applies GQ and DP as per-variant medians across samples, not per-genotype
# cutoffs. Completeness >= 0.9 is expressed as AN >= 0.9 * 2N.
N_SAMPLES=$(bcftools query -l "${OUT_PREFIX}_s1_biallelic.vcf.gz" | wc -l)
MAX_AN=$((N_SAMPLES * 2))
log "Samples: $N_SAMPLES"

bcftools filter \
    --exclude \
    "MEDIAN(FMT/GQ) < ${MIN_GQ} || \
     MEDIAN(FMT/DP) < ${MIN_DP} || \
     AN < (${MAX_AN} * 0.9)" \
    "${OUT_PREFIX}_s1_biallelic.vcf.gz" \
    -O z -o "${OUT_PREFIX}_s2_median_filters.vcf.gz"

bcftools index --tbi "${OUT_PREFIX}_s2_median_filters.vcf.gz"
log "Variants after median GQ/DP and completeness: $(bcftools view -H ${OUT_PREFIX}_s2_median_filters.vcf.gz | wc -l)"

# AB ratio = fraction of het calls whose REF read fraction is within [0.15, 0.85].
# Sites with at least 5 hets and AB ratio < 0.9 fail. Sites without AD are kept.
python3 - \
    "${OUT_PREFIX}_s2_median_filters.vcf.gz" \
    "${OUT_PREFIX}_ab_fail_sites.txt" \
    <<'PYEOF'
import sys, gzip

vcf_in  = sys.argv[1]
out_txt = sys.argv[2]

AB_THRESHOLD = 0.9
BALANCE_MIN  = 0.15
BALANCE_MAX  = 0.85

fail_sites = []
opener = gzip.open if vcf_in.endswith('.gz') else open

with opener(vcf_in, 'rt') as fh:
    for line in fh:
        if line.startswith('#'):
            continue
        fields  = line.rstrip('\n').split('\t')
        chrom, pos, vid = fields[0], fields[1], fields[2]
        fmt     = fields[8].split(':')
        samples = fields[9:]

        try:
            gt_idx = fmt.index('GT')
            ad_idx = fmt.index('AD')
        except ValueError:
            continue

        n_het = 0
        n_imbalanced = 0
        for samp in samples:
            vals = samp.split(':')
            if len(vals) <= max(gt_idx, ad_idx):
                continue
            gt = vals[gt_idx].replace('|', '/')
            if gt in ('./.', '.'):
                continue
            alleles = gt.split('/')
            if len(alleles) == 2 and alleles[0] != alleles[1] and '.' not in alleles:
                n_het += 1
                ad_str = vals[ad_idx]
                if ',' not in ad_str or '.' in ad_str:
                    continue
                try:
                    ad = [int(x) for x in ad_str.split(',')]
                    total = sum(ad)
                    if total == 0:
                        continue
                    ref_frac = ad[0] / total
                    if ref_frac < BALANCE_MIN or ref_frac > BALANCE_MAX:
                        n_imbalanced += 1
                except ValueError:
                    continue

        if n_het >= 5:
            ab_ratio = 1.0 - (n_imbalanced / n_het)
            if ab_ratio < AB_THRESHOLD:
                fail_sites.append(f"{chrom}\t{pos}\t{vid}")

print(f"Sites failing AB ratio filter: {len(fail_sites)}")
with open(out_txt, 'w') as out:
    out.write('\n'.join(f.split('\t')[2] for f in fail_sites) + '\n')
PYEOF

# The fail list is matched on the ID column, so the VCF needs unique variant IDs.
N_AB_FAIL=$(wc -l < "${OUT_PREFIX}_ab_fail_sites.txt" || echo 0)
if [[ "$N_AB_FAIL" -gt 0 ]]; then
    bcftools view \
        --exclude "ID=@${OUT_PREFIX}_ab_fail_sites.txt" \
        "${OUT_PREFIX}_s2_median_filters.vcf.gz" \
        -O z -o "${OUT_PREFIX}_s3_ab_filtered.vcf.gz"
else
    cp "${OUT_PREFIX}_s2_median_filters.vcf.gz" "${OUT_PREFIX}_s3_ab_filtered.vcf.gz"
fi

bcftools index --tbi "${OUT_PREFIX}_s3_ab_filtered.vcf.gz"
log "Variants after AB filter: $(bcftools view -H ${OUT_PREFIX}_s3_ab_filtered.vcf.gz | wc -l)"

# MAF here is in the study data only; MAF in 1KGP3 is applied in step 3.
plink \
    --vcf "${OUT_PREFIX}_s3_ab_filtered.vcf.gz" \
    --make-bed \
    --snps-only just-acgt \
    --biallelic-only strict \
    --maf "$MAF" \
    --geno "$MAX_MISS" \
    --autosome \
    --keep-allele-order \
    --threads "$THREADS" \
    --out "${OUT_PREFIX}_s4_plink_qc"

log "Variants after MAF/missingness: $(wc -l < ${OUT_PREFIX}_s4_plink_qc.bim)"

awk '($5=="A" && $6=="T") || ($5=="T" && $6=="A") || \
     ($5=="C" && $6=="G") || ($5=="G" && $6=="C") {print $2}' \
    "${OUT_PREFIX}_s4_plink_qc.bim" > "${OUT_PREFIX}_atgc_snps.txt"

log "AT/GC SNPs removed: $(wc -l < ${OUT_PREFIX}_atgc_snps.txt)"

plink \
    --bfile "${OUT_PREFIX}_s4_plink_qc" \
    --exclude "${OUT_PREFIX}_atgc_snps.txt" \
    --make-bed \
    --keep-allele-order \
    --threads "$THREADS" \
    --out "${OUT_PREFIX}_noatgc"

log "Final variant count: $(wc -l < ${OUT_PREFIX}_noatgc.bim) (${OUT_PREFIX}_noatgc)"
