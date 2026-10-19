#!/usr/bin/env bash
# Merge the cases with the HC, DEE and GEL Epi control sets (stepwise, on shared SNPs).
# Usage: export PIPELINE_CONFIG=/path/to/config.sh; bsub < 01_merge_controls.sh
#BSUB -P re_gecip_neurology
#BSUB -q medium -M 16000 -n 2 -R "rusage[mem=16000]"
#BSUB -J prs01 -o logs/%J_prs01.out -e logs/%J_prs01.err
set -euo pipefail
source "${PIPELINE_CONFIG:?export PIPELINE_CONFIG=/path/to/config.sh}"

cd ${WORKDIR}
mkdir -p merged
CASE=${CASE_PLINK}
WORK=${WORKDIR}/merged/_merge_tmp
mkdir -p ${WORK}

awk '{print $1"\t"$2}' ${CASE}.fam > merged/case_samples.txt

merge_in() {
    local target=$1 add=$2 label=$3 out=$4
    echo "  merging ${label}"

    # a case that is also in this control cohort stays a case only
    awk 'NR==FNR{c[$2]; next} ($2 in c){print $1"\t"$2}' \
        ${CASE}.fam ${add}.fam > merged/${label}_overlap_with_cases.txt
    local n_overlap=$(wc -l < merged/${label}_overlap_with_cases.txt | tr -d ' ')
    local remove=""
    if [[ $n_overlap -gt 0 ]]; then
        remove="--remove merged/${label}_overlap_with_cases.txt"
    fi
    echo "  ${label}: ${n_overlap} cases also in this cohort, dropped from the controls"

    awk '{print $2}' ${target}.bim | sort -u > ${WORK}/t.txt
    awk '{print $2}' ${add}.bim    | sort -u > ${WORK}/a.txt
    comm -12 ${WORK}/t.txt ${WORK}/a.txt > ${WORK}/keep.txt

    module load ${PLINK2_MODULE}
    plink2 --bfile ${target} --extract ${WORK}/keep.txt --make-bed --out ${WORK}/t_sub
    plink2 --bfile ${add}    --extract ${WORK}/keep.txt ${remove} --make-bed --out ${WORK}/a_sub

    module unload ${PLINK2_MODULE}
    module load ${PLINK19_MODULE}

    set +e
    plink --bfile ${WORK}/t_sub \
          --bmerge ${WORK}/a_sub.bed ${WORK}/a_sub.bim ${WORK}/a_sub.fam \
          --make-bed --out ${out}
    set -e

    # strand-flip retry
    if [ -s ${out}-merge.missnp ]; then
        plink --bfile ${WORK}/a_sub --flip ${out}-merge.missnp \
              --make-bed --out ${WORK}/a_flip
        plink --bfile ${WORK}/t_sub \
              --bmerge ${WORK}/a_flip.bed ${WORK}/a_flip.bim ${WORK}/a_flip.fam \
              --make-bed --out ${out}
    fi
    module unload ${PLINK19_MODULE}

    awk '{print $1"\t"$2}' ${WORK}/a_sub.fam > merged/${label}_samples.txt
    echo "  ${label} n=$(wc -l < merged/${label}_samples.txt), merged fam=$(wc -l < ${out}.fam)"
}

merge_in ${CASE}            ${HC_PLINK}  hc  merged/case-HC
merge_in merged/case-HC     ${DEE_PLINK} dee merged/case-HC-DEE
merge_in merged/case-HC-DEE ${EPI_PLINK} epi merged/case-HC-DEE-EPI

echo "Final merged: $(wc -l < merged/case-HC-DEE-EPI.fam) samples, $(wc -l < merged/case-HC-DEE-EPI.bim) variants"
