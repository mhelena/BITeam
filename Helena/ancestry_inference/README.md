# Ancestry inference

Snakemake pipeline that assigns a continental super-population (AFR, AMR, EAS, EUR, SAS) to
each sample in a multi-sample VCF, following the Genomics England AggV2 ancestry protocol
(https://re-docs.genomicsengland.co.uk/ancestry_inference/).

Steps:

0. Optional: merge per-sample VCFs (`scripts/00_merge_vcfs.sh`).
1. Variant QC of the study VCF: autosomal biallelic SNPs, per-variant median GQ and DP,
   completeness, allelic balance, MAF, missingness, AT/GC removal (`01_qc_vcf.sh`).
2. Download 1000 Genomes phase 3 (1KGP3) and write sample and super-population lists
   (`02_prepare_reference.sh`).
3. Merge with 1KGP3, harmonise alleles, drop high-LD regions, LD prune on reference samples,
   HWE within each of AFR/EAS/EUR/SAS (`03_merge_and_prune.sh`).
4. KING check on the reference, PCA on unrelated 1KGP3, projection of study samples
   (`04_pca_gcta.sh`).
5. Random forest on PCs 1-8 of the reference; study samples are assigned when the top
   probability is >= 0.8, with a kNN fallback for the rest (`05_ancestry_classifier.py`).
6. Plots (`06_plot_results.py`).

## Requirements

```
conda env create -f environment.yml
conda activate ancestry_inference
```

`bin/` must be on `PATH` (run_pipeline.sh does this). It holds `gcta64`, which is not GCTA:
it is a small Python wrapper around plink2 that mimics `gcta64 --pca` and
`gcta64 --project-loading`, kept under that name so step 4 calls it unchanged. It uses the
`plink2` found on `PATH`, or the one given in the `PLINK2` environment variable.

## Configure

`config/config.yaml` is the local profile and holds every setting. `config/cluster.yaml` is
layered on top of it for HPC runs and only lists what differs (input as a list of per-sample
VCFs, logs under `results/`, median GQ and DP of 30).

Input, set one of:

- `input_vcf`: one multi-sample VCF
- `input_vcf_dir`: a directory of per-sample VCFs
- `input_vcf_list`: a list of per-sample VCF paths

Other keys: `genome_build` (hg38 or hg19), `ref_1kg_plink` (prefix of an existing 1KGP3
PLINK set, see notes), `output_dir`, `log_dir`, `qc`, `ld_pruning`, `pca`, `relatedness`,
`classifier`, `threads`, `memory`.

## Run

From this directory:

```
bash run_pipeline.sh local            # or: bash run_pipeline.sh local -n  for a dry run
mkdir -p logs && qsub run_pipeline.sh cluster
```

On SGE, set `CONDA_SH` (path to conda.sh) and `CONDA_ENV` so the job can activate the
environment, e.g. `qsub -v CONDA_SH=...,CONDA_ENV=ancestry_inference run_pipeline.sh cluster`.
Extra arguments after the profile go straight to snakemake.

## Outputs

Under `output_dir`:

- `05_classify/ancestry_ancestry_assignments.tsv`: one row per study sample with RF
  probabilities, PCs, `MaxProb`, `PredLabel`, the final `AncestryGEL` label and
  `ClassificationMethod` (RF, kNN or Unassigned)
- `05_classify/ancestry_summary.json`: OOB error, CV accuracy, counts per group
- `05_classify/ancestry_reference_predictions.tsv`, `ancestry_rf_model.joblib`
- `04_pca/_all_pcs.tsv`: reference and projected study PCs
- `06_plots/`: PC1/PC2, PC1/PC3 and PC2/PC3 scatters, probability boxplots, counts, and the
  confusion matrix on reference samples

## Notes

- PCA is computed on the reference only and study samples are projected. The wrapper
  re-scores the reference samples with `plink2 --score`, the same way the study samples are
  projected, so both sets of PCs are on one scale.
- Samples below the RF threshold go to a kNN vote among the 50 nearest reference samples and
  are assigned if at least 70% agree. This rescues samples pushed off the reference clusters
  by batch effects.
- The 5-fold CV in step 5 runs with `n_jobs=1`; in parallel it ran out of memory.
- `ld_pruning.r2` is 0.5 and the local profile uses median GQ 15 and DP 10. GEL uses
  r2 = 0.1 and 30/30. These are the values the pipeline has been run with.
- The allelic balance filter matches failing sites by the VCF ID column, so variant IDs
  should be unique (not `.`).
- Step 2 downloads 1KGP3 but its output is not directly usable by step 3 (the hg38 `.pvar`
  stays zst-compressed and the `.psam` has a different name; the hg19 VCFs are never
  converted). In practice set `ref_1kg_plink` to a prepared PLINK set. Step 2 still runs
  (and downloads) because steps 3 to 6 use its sample and label lists.
- Step 0 merges with `--missing-to-ref`, so a site absent from one sample's VCF becomes 0/0.

## Extra: genetic similarity

`extras/genetic_similarity.R` estimates proportional similarity to 21 reference populations
and an approximate FST to each (Prive et al. 2022, bigsnpr). It runs on any PLINK bed set,
for example the step 1 output `01_qc/study_qc_noatgc`:

```
Rscript extras/genetic_similarity.R results/01_qc/study_qc_noatgc /path/to/ref_dir results/gen_sim
```

`ref_dir` must contain `ref_freqs.csv.gz` and `projection.csv.gz` from the bigsnpr ancestry
tutorial (https://privefl.github.io/bigsnpr/articles/ancestry.html). Needs R >= 4.1 with
bigsnpr, bigreadr, quadprog and Matrix.
