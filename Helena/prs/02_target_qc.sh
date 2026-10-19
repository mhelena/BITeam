#!/usr/bin/env bash
# Joint QC of the merged set: MAF, HWE, missingness, heterozygosity, relatedness.
# Usage: bsub < 02_target_qc.sh
#BSUB -P re_gecip_neurology
#BSUB -q medium -M 16000 -n 2 -R "rusage[mem=16000]"
#BSUB -J prs02 -o logs/%J_prs02.out -e logs/%J_prs02.err
set -euo pipefail
source "${PIPELINE_CONFIG:?export PIPELINE_CONFIG=/path/to/config.sh}"

module load ${PLINK2_MODULE}
cd ${WORKDIR}

IN=merged/case-HC-DEE-EPI
OUT=merged/case-HC-DEE-EPI.QC
echo "input: $(wc -l < ${IN}.fam) samples, $(wc -l < ${IN}.bim) variants"

plink2 --bfile ${IN} \
       --maf 0.01 --hwe 1e-6 --geno 0.01 --mind 0.01 \
       --write-snplist --make-just-fam \
       --threads 2 --memory 14000 --out ${OUT}
echo "MAF/HWE/geno/mind: $(wc -l < ${OUT}.fam) samples, $(wc -l < ${OUT}.snplist) variants"

# pruned set is only used for the heterozygosity and relatedness checks
plink2 --bfile ${IN} --keep ${OUT}.fam --extract ${OUT}.snplist \
       --indep-pairwise 200 50 0.25 \
       --threads 2 --memory 14000 --out ${OUT}

plink2 --bfile ${IN} --keep ${OUT}.fam --extract ${OUT}.prune.in \
       --het --threads 2 --memory 14000 --out ${OUT}

awk 'NR==1 {for (i=1; i<=NF; i++) {if ($i=="F") f=i; if ($i=="IID") id=i}; next}
     {F[$id]=$f; s+=$f; ss+=$f*$f; n++}
     END {m=s/n; sd=sqrt(ss/n-m*m)
          for (k in F) if (F[k] < m-3*sd || F[k] > m+3*sd) print k}' \
    ${OUT}.het > ${OUT}.het_outliers.txt
awk 'NR==FNR {out[$1]; next} !($2 in out) {print $1"\t"$2}' \
    ${OUT}.het_outliers.txt ${OUT}.fam > ${OUT}.valid_after_het.txt
echo "heterozygosity (mean +/- 3 SD): $(wc -l < ${OUT}.het_outliers.txt | tr -d ' ') removed"

# KING kinship > 0.0884 is 3rd degree or closer
plink2 --bfile ${IN} --keep ${OUT}.valid_after_het.txt --extract ${OUT}.prune.in \
       --king-cutoff 0.0884 \
       --threads 2 --memory 14000 --out ${OUT}.king
KEEP=${OUT}.king.king.cutoff.in.id
[[ -s ${KEEP} ]] || KEEP=${OUT}.king.cutoff.in.id
n_rel=$(( $(wc -l < ${OUT}.valid_after_het.txt) - $(awk 'NR>1' ${KEEP} | wc -l) ))
echo "relatedness: ${n_rel} removed"

plink2 --bfile ${IN} --keep ${KEEP} --extract ${OUT}.snplist \
       --make-bed --threads 2 --memory 14000 --out ${OUT}

echo "final: $(wc -l < ${OUT}.fam) samples, $(wc -l < ${OUT}.bim) variants"
