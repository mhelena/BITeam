# Case-control PRS in the GEL Research Environment

PRSice-2 polygenic scores for a case cohort against GEL control cohorts (healthy
controls, DEE cases, epilepsy controls), on LSF in the Research Environment. The
input is a QC'd PLINK set of cases, as produced by `gel_cohort_extraction`, plus the
three control sets.

## Setup

1. Copy `config.sh`, set `PROJECT_DIR`, `COHORT` and the input paths.
2. `export PIPELINE_CONFIG=/absolute/path/to/config.sh` in the shell you submit from.
   LSF passes the environment to the job, which is how the scripts find the config
   (`bsub < script.sh` runs a spooled copy, so the script cannot locate files next to itself).
3. `mkdir -p logs` in the folder you submit from.
4. The `#BSUB -P` lines use `re_gecip_neurology`; change them if your project differs.

R scripts take `WORKDIR` as the first argument or read the settings from the
environment, so run them as `source $PIPELINE_CONFIG && Rscript script.R`, or with
`bsub -P $LSF_PROJECT -q short -M 8000 "module load $R_MODULE && Rscript script.R"`.
In RStudio, edit the default paths at the top of the script and Source it.

## Core, in order

| Script | What it does | Run |
|---|---|---|
| `01_merge_controls.sh` | Merge cases with HC, DEE and GEL Epi on shared SNPs; cases found in a control set are dropped from that set | `bsub <` |
| `02_target_qc.sh` | Joint QC: MAF 0.01, HWE 1e-6, geno/mind 0.01, heterozygosity mean +/- 3 SD, KING 0.0884 | `bsub <` |
| `03_pca_vs_1kg.sh`, `03b_plot_pca_1kg.R` | Optional ancestry check against 1000 Genomes | `bsub <`, then Rscript |
| `04_pca_target_only.sh` | 10 PCs of the QC'd target | `bsub <` |
| `05_build_cov_pheno.sh` | `cov.txt` (IID PC1-4 Sex) and `pheno.txt` (cases 1, `PRSICE_CONTROL` 0) | `bash` |
| `06_run_prsice.sh` | PRSice for each score in `PRS_LIST`; every sample is scored | `bsub <` |
| `07_plot_boxplots.R` | Boxplots by cohort, case vs each control from `lm(PRS_z ~ cohort + PC1-4 + Sex)` | Rscript |

PRSice runs with `--ignore-fid`, so the covariate and phenotype files are keyed on
IID. Only the cases and `PRSICE_CONTROL` are in `pheno.txt`, so the other two control
cohorts do not influence the threshold PRSice picks. Downstream scripts take each
score at the best threshold from `PRS_<name>.summary`.

Outputs, all under `WORKDIR`: `merged/case-HC-DEE-EPI.QC.*` and the per-cohort sample
lists, `pca/`, `cov.txt`, `pheno.txt`, `results/PRS_<name>.*` and the boxplots.

## Modules

All run after the core and write next to its outputs.

- `modules/cohort_factor_model.R`: one model per score, `PRS_z ~ Cohort + PC1-4 + Sex`,
  cases vs each control; table and forest plot in `results/`.
- `modules/summarize_prsice_thresholds.R [results_dir] [suffix]`: Num_SNP and R2 at
  every threshold, best Pt marked. Works on any PRSice results folder.
- `modules/apoe/`: APOE sensitivity for the longevity score (`LONGEVITY_SCORE`, which
  must be in `PRS_LIST`). `01` calls APOE genotypes, `02` tests the score by APOE
  genotype and cohort, `03` reruns PRSice with chr19:44.4-45.5 Mb removed from the base
  GWAS, `04` compares the cohort P values with and without APOE. Run `02` after `01`,
  and `04` after `03`.
- `modules/within_case/`: scores tested against a binary outcome measured in the cases
  (`OUTCOME`, `OUTCOME_TSV`), and optionally against age at onset (`ONSET_TSV`).
  `01` subsets and re-filters the cases, `02` computes PCs within them, `03` writes the
  phenotype and covariate files, `04` runs PRSice, `05` gives the summary table and
  plots, `06` fits all `JOINT_PRS_LIST` scores in one logistic model, and `07` repeats
  that model with thresholds fixed from the case-control run, repeated 5-fold CV with
  shared folds, label permutation, PRS-block residual permutation and, if
  `EXTRA_COV_TSV` exists, the same tests with those covariates in the baseline.
  `EXTRA_COV_TSV` is tab-separated with IID plus any columns; numbers are treated as
  measurements, so give categories (centre, batch) as text.

R packages: ggplot2, dplyr, patchwork (APOE comparison only).
