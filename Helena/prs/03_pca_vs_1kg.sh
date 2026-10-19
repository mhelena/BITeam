#!/usr/bin/env bash
# Optional ancestry check: PCA of the QC'd target merged with 1000 Genomes.
# Usage: bsub < 03_pca_vs_1kg.sh, then Rscript 03b_plot_pca_1kg.R
#BSUB -P re_gecip_neurology
#BSUB -q short -M 12000 -n 2 -R "rusage[mem=12000]"
#BSUB -J prs03 -o logs/%J_prs03.out -e logs/%J_prs03.err
set -euo pipefail
source "${PIPELINE_CONFIG:?export PIPELINE_CONFIG=/path/to/config.sh}"

cd ${WORKDIR}
mkdir -p pca
TARGET=merged/case-HC-DEE-EPI.QC
WORK=${WORKDIR}/pca/_tmp
mkdir -p ${WORK}

# 1KG SNP IDs use dashes, ours use colons
awk 'BEGIN{OFS="\t"} {gsub("-",":",$2); print}' ${KG_PLINK}.bim > ${WORK}/kg.bim
cp ${KG_PLINK}.bed ${WORK}/kg.bed
cp ${KG_PLINK}.fam ${WORK}/kg.fam

module load ${PLINK2_MODULE}

awk '{print $2}' ${TARGET}.bim | sort -u > ${WORK}/t.txt
awk '{print $2}' ${WORK}/kg.bim | sort -u > ${WORK}/kg.txt
comm -12 ${WORK}/t.txt ${WORK}/kg.txt > ${WORK}/keep.txt

plink2 --bfile ${TARGET}  --extract ${WORK}/keep.txt --make-bed --out ${WORK}/t_sub
plink2 --bfile ${WORK}/kg --extract ${WORK}/keep.txt --make-bed --out ${WORK}/kg_sub

module unload ${PLINK2_MODULE}
module load ${PLINK19_MODULE}

set +e
plink --bfile ${WORK}/t_sub \
      --bmerge ${WORK}/kg_sub.bed ${WORK}/kg_sub.bim ${WORK}/kg_sub.fam \
      --make-bed --out pca/target-KG
set -e
if [ -s pca/target-KG-merge.missnp ]; then
    plink --bfile ${WORK}/kg_sub --flip pca/target-KG-merge.missnp \
          --make-bed --out ${WORK}/kg_flip
    plink --bfile ${WORK}/t_sub \
          --bmerge ${WORK}/kg_flip.bed ${WORK}/kg_flip.bim ${WORK}/kg_flip.fam \
          --make-bed --out pca/target-KG
fi

module unload ${PLINK19_MODULE}
module load ${PLINK2_MODULE}

plink2 --bfile pca/target-KG --indep-pairwise 200 50 0.1 --out pca/target-KG_prune
plink2 --bfile pca/target-KG \
       --extract pca/target-KG_prune.prune.in \
       --pca 20 --out pca/pca_target_KG

echo "Wrote pca/pca_target_KG.{eigenvec,eigenval}"
