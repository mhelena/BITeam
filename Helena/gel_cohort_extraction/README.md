# GEL cohort extraction

Pulls a set of Genomics England participants out of the aggV2 aggregate VCF chunks
into a QC'd PLINK BED, then runs per-cohort sample QC. Runs inside the GEL Research
Environment and submits jobs with LSF (`bsub`).

## What it does

1. `01_eur_filter.sh` keeps samples with predicted European ancestry >= 0.8 and
   converts participant IDs to platekeys.
2. `02_extract_plink/` runs one job per aggV2 chunk: FILTER=PASS, biallelic SNPs,
   autosomes, `medianDepthNonMiss >= 10`, `medianGQ >= 15`, subset to the cohort,
   missing genotypes set to reference, `AC > 0`, then per-site allele balance
   (binomial p < 0.01, at least 25% of hets balanced), missingness (<= 5% DP=0) and
   completeness (>= 50% called) in `_ab_miss_complete.py`. Chunks are merged per
   chromosome and then genome-wide, where MAF >= 0.05, geno < 0.01 and removal of
   A/T and C/G SNPs are applied.
3. `03_sample_qc/` removes samples with Ti/Tv outside mean +/- 3 SD and variants
   failing HWE at p < 1e-5.

`build_cohorts.R` is optional. It builds four control cohorts from LabKey
(epilepsy controls, DEE, healthy candidates, SCN1A screening pool), all EUR >= 0.8
and GRCh38, and writes a `samples.txt` for each.

## Configure

Copy `config.sh` somewhere on the RE, set `PROJECT_DIR` (and `LSF_PROJECT` if you
are not in the neurology GeCIP), then:

```bash
export PIPELINE_CONFIG=/absolute/path/to/config.sh
```

Every script reads the config through `PIPELINE_CONFIG` and stops with an error if it
is unset. An environment variable is used because `bsub < script.sh` runs a spooled
copy of the script, so `$(dirname $0)` does not point at this folder. The submit
scripts pass the variable on to the jobs with `-env "all"`.

`COHORT` is read from the environment if set, so one config can serve several
cohorts (`COHORT=dee ./02_extract_plink/submit.sh`).

Two reference files go under `${PROJECT_DIR}/reference/`:

- `chunks.txt`: one aggV2 chunk name per line, for example
  `ls ${AGGV2_DIR}/gel_mainProgramme_aggV2_*.vcf.gz | sed 's/.*aggV2_//; s/\.vcf\.gz$//'`
- `aggV2_ancestry_sex.tsv`: tab-separated with a header, columns `platekey`,
  `participant_id`, `pred_european_ancestries`, `sex` (from the LabKey table
  `aggregate_gvcf_sample_stats`).

## Run order

```bash
# optional, writes ${PROJECT_DIR}/<cohort>/samples.txt for the control cohorts
module load R/4.3.3
Rscript build_cohorts.R "$PROJECT_DIR" [non_epi_terms.txt]

# otherwise put your own list in ${PROJECT_DIR}/${COHORT}/samples.txt
bash 01_eur_filter.sh
./02_extract_plink/submit.sh          # optional queue argument, default long
./03_sample_qc/submit.sh              # after the stage 2 jobs finish (bjobs)
```

`build_cohorts.R` uses `LABKEY_FOLDER` (default: main programme v19) and
`LABKEY_BASE` if set. The optional term file lists disease or HPO terms, one per line
(Supplementary Table 3 of PMID:37006128); healthy candidates matching any of them
are dropped.

Stage 2 jobs skip chunks that already have output, so resubmitting after a failure
only reruns the missing ones.

## Outputs

Under `${PROJECT_DIR}/${COHORT}/`:

- `samples_eur08.txt`: platekeys kept by stage 1
- `final/${COHORT}_ALL.{bed,bim,fam}`: extracted and variant-QC'd cohort
- `qc/${COHORT}_qc.{bed,bim,fam}`: after Ti/Tv and HWE, with `qc/${COHORT}_titv_summary.txt`,
  `qc/${COHORT}_qc_summary.txt` and the HWE report `qc/${COHORT}_qc.hardy`
- `logs/`: LSF logs
- `scratch/`: per-chunk and per-chromosome intermediates, safe to delete once
  `final/` exists
