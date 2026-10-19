# Epilepsy genomics pipelines

Analysis pipelines used by the group, mostly for whole-genome data in the Genomics England
Research Environment (GEL RE, LSF), the UCL CS cluster (SGE) or a local machine.

| Folder | What it does | Runs on |
|---|---|---|
| `ancestry_inference/` | Genetic ancestry from a multi-sample VCF: QC, merge with 1000 Genomes, PCA, random forest classifier | local or HPC (Snakemake) |
| `gel_cohort_extraction/` | Pull a cohort out of the GEL aggV2 VCF chunks into a QC'd PLINK set, then sample QC | GEL RE |
| `gene_carrier_screen/` | Find participants carrying unique missense variants in a gene of interest, from disease categories unrelated to the phenotype | GEL RE |
| `panel_carrier_exclusion/` | Remove carriers of rare variants in a PanelApp gene panel from a control pool | GEL RE |
| `prs/` | Case-control PRS with PRSice-2: merge, target QC, PCA, covariates, scores, plots, contrasts, plus optional modules | GEL RE |
| `rare_variant_vcf_qc/` | Variant QC of a joint-called VCF for rare-variant work, and building a multi-cohort VCF for burden tests | UCL CS, local, GEL RE |
| `brainage_prs/` | PRS against predicted brain-age gap in ENIGMA-Epilepsy | local |
| `ldsc_genetic_correlation/` | Pairwise genetic correlations between GWAS summary statistics with LDSC, plus heatmap and forest plot | local or HPC |

Each folder has its own README with inputs, configuration and run order.

## Data

No data are kept in this repository. Individual-level data, sample lists, phenotype and
covariate files stay where they were generated (for GEL data, inside the Research
Environment). `.gitignore` blocks the usual genotype, score and spreadsheet formats; check
`git status` before every commit anyway.

## Conventions

- Paths and project settings live in each pipeline's config file, never in the scripts.
- LSF scripts are submitted with `bsub < script.sh`. LSF runs a spooled copy, so scripts find
  their config through an absolute path or an environment variable, not `$(dirname $0)`.
- Keep comments short and only where the reason is not obvious from the code.
