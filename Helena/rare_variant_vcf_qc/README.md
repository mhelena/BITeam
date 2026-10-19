# Rare-variant VCF QC

Site-level QC of a joint-called, PASS-only, GRCh38 multi-sample VCF for rare-variant work,
and a small LSF pipeline that turns per-cohort PLINK sets from GEL into one merged VCF for
burden tests.

```
helpers/      ab_miss_complete.py, hwe_exact.py (used by both QC runners)
sge/          QC on the UCL CS cluster (qsub), one array job per chromosome
local/        the all-ancestry QC on a single machine
gel_burden/   multi-cohort burden VCF on the GEL Research Environment (bsub)
```

Each runner reads the `config.sh` in its own folder, or the file named in `PIPELINE_CONFIG`.

## Filters

A site is removed if any of these hold:

- fewer than 25% of heterozygous calls are allele-balanced (binomial p >= 0.01 around 0.5)
- more than 5% of samples have DP = 0
- fewer than 50% of samples have a called genotype
- the indel is longer than 50 bp (`abs(strlen(REF)-strlen(ALT))>50`)

plus, depending on the mode:

| | `allanc` | `eur` |
|---|---|---|
| Samples | all | IIDs in `EUR_KEEP` |
| MEDIAN_DP < 10 or MEDIAN_GQ < 15 | yes | no |
| HWE exact p < 1e-5 | no | yes |
| Chromosomes | every chr1-22, X, Y, M contig in the header | same |

HWE is left out of the all-ancestry mode because pooling ancestries inflates HWE deviation
(Wahlund effect) and it adds little for rare variants. MEDIAN_DP and MEDIAN_GQ have to be
in INFO already; the scripts do not compute them.

There is no MAF, LD or allele-count filter, and multi-allelic sites are kept.

The repeat/low-complexity (LCR) mask is skipped on purpose in the all-ancestry output, so
users can apply their own. On SGE, setting `LCR_BED` in `sge/config.sh` (for example GIAB
`GRCh38_AllTandemRepeatsandHomopolymers_slop5.bed.gz`, tabix-indexed) adds it as an extra
exclusion in the filter step for either mode.

## SGE (UCL CS)

```
cd sge
# edit config.sh: WORKDIR, RAW_VCF, EUR_KEEP (eur only), tool paths
./qc_submit.sh                  # MODE from config.sh, allanc by default
MODE=eur ./qc_submit.sh
```

Stages: 0 split by chromosome, 1 fail sites per chromosome, 2 union of fails,
3 filter per chromosome, 4 concatenate. Output in `OUTDIR`:
`study_qc_<mode>.vcf.gz` (+ `.tbi`), `study_qc_<mode>.failsites.tsv` (CHROM, POS, reason)
and `study_qc_<mode>.filters.log`.

Per-chromosome stages write `.done` flags, so running `qc_submit.sh` again after a failure
only redoes what is missing. `./qc_submit.sh --redo chr10,chr16` also deletes everything
already written for those chromosomes first. `OUTDIR/chrom_list.txt` is built on the first
run; edit it to process a subset.

bcftools 1.9 is enough. Its `+fill-tags` has no HWE, which is why `hwe_exact.py` exists.

## Single machine

```
cd local
# edit config.sh: IN, OUTDIR, OUT_PREFIX
./qc_local_run.sh               # both passes, several hours on a 50-60 GB VCF
./qc_local_run.sh --pass2-only  # rerun the filter using pass 1's all_fail_regions.tsv.gz
```

This is the all-ancestry mode only. Pass 1 computes the allele balance, missingness and
completeness fails; pass 2 removes them and applies MEDIAN_DP/GQ and the indel limit as a
`bcftools view -e` expression. Unlike the SGE version, a site with MEDIAN_DP or MEDIAN_GQ
missing (".") is kept here, since a missing value never satisfies `INFO/MEDIAN_DP<10` in
bcftools.

## GEL RE burden VCF

Input is one PLINK set per chromosome per cohort, as written by `gel_cohort_extraction`
(`<cohort>_chr<N>.{bed,bim,fam}`).

```
cd gel_burden
# edit config.sh: PROJECT_DIR, COHORTS, MERGE_OUT_TAG, perchr_dir()
./burden_submit.sh build        # per-cohort VCFs, then bcftools merge -m none
# run variant QC on ${MERGE_OUT_TAG}_for_burden.vcf.gz, giving ACGT0_IN
./burden_submit.sh acgt0        # recompute AC/AN/AF, drop AC=0 sites
```

The variant QC between `build` and `acgt0` is not part of these scripts. Set
`OUTPUT_COPY_DIR` to have the final VCF copied to a shared folder and checked by md5 and
record count.

On GEL RE, `bcftools/1.16` does not include bgzip or tabix; they come from the `htslib/1.21`
module (not samtools). Indexes are written with `bcftools index --tbi`.
