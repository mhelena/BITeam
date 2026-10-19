# PRS and brain-age gap in ENIGMA-Epilepsy

Tests whether polygenic risk scores are associated with predicted brain-age gap (pBAG, also
called BrainPAD) in people with epilepsy. PRSice-2 scores 8 GWAS bases against 4 clinical
contrasts, then the selected PRS are regressed on pBAG alongside clinical covariates, and all
tests are corrected with Benjamini-Hochberg FDR.

Cases whose MRI finding is "Post surgery" or "Large lesion" are excluded in every step.

## Steps

| Step | Script | What it does |
|---|---|---|
| 0 | `00_create_pheno_files.R` | Binary phenotype files for the 4 PRSice contrasts |
| 1 | `01_M1_run_PRSice.sh` | PRSice-2, 8 GWAS x 4 contrasts (32 runs) |
| 2 | `02_M1_summarise_PRSice.R` | One table of the 32 best-threshold results |
| 3 | `03_M2_clinical_only.R` | Clinical-only baseline model of pBAG |
| 4 | `04_M3_individual_PRS.R` | pBAG ~ PRS + PC1-4 + Sex for 5 PRS x 7 thresholds; picks the best threshold per PRS |
| 5 | `05_M4_hierarchical.R` | Block regression: PCs, + demographics, + epilepsy factors, + 5 PRS |
| 6 | `06_M5_stratification.R` | PRS tertiles/quartiles (ANOVA, Kruskal-Wallis, trend) and PRS x clinical interactions |
| 7 | `07_FDR_correction.py` | BH-FDR within each test family and across all 85 tests |

The 5 PRS used from step 4 on are Epilepsy, Intelligence, Longevity, Frailty and Brain Age
(GWM), chosen from the 8 bases by LDSC genetic correlation. Steps 4-6 drop or set to NA any
PRS more than 3 SD from the mean. Steps 5 and 6 read the best thresholds from the step 4
output.

Step 7 reads the M3 p-values from the step 4 output, but the M1, M4 and M5 p-values are
typed into the script from the reference run. Update them by hand if you re-run steps 2, 5
or 6 on new data.

## Inputs

All files are tab-separated with a header.

Clinical metadata (`METADATA`), one row per participant:
- `FID`, `IID`: PLINK IDs
- `MRI_findings`: rows with "Post surgery" or "Large lesion" are dropped
- `pBAG_winsorized`: predicted brain-age gap, winsorised at +/- 2 SD
- `AgeAtFirstScan`: age at first MRI (years)
- `Age at onset`: age at first seizure (years)
- `Years.Of.Epilepsy`: duration; equal to age at scan minus age at onset, so not modelled
- `EpilepsySyndrome_forPRS`: "Focal" / "Generalised" (other values are treated as missing)
- `DRE`: "DRE at time of scan" / "NO DRE at time of scan"
- `HistoryOfConvulsions`: "History of convulsions" / "NO History of convulsions"
- `VPA`: "VPA exposure" / "NO VPA exposure"

Covariates (`COVARIATES`): `FID`, `IID`, `PC1`-`PC4`, `Sex`. The PCs are the first 4
ancestry components from a PLINK PCA of the QC'd European target set.

PLINK target genotypes (`TARGET`): the `.fam` is also used in step 0 to map `FID` to the
genotype IDs.

GWAS base files (8, GRCh37) with columns `SNP`, `CHR`, `BP`, `A1`, `A2`, `BETA`, `P`: ILAE3
epilepsy, intelligence (Savage 2018), longevity (Deelen 2019), Alzheimer's disease (Kunkle
2019), ENIGMA brain age (GM, WM, GWM) and frailty (Atkins 2021).

`VALID_SNPS`: the `.valid` SNP list from an earlier PRSice run, passed to `--extract` so all
runs use the same SNP set.

Files written by step 0 (`PHENO_DIR`), each with `FID`, `IID`, `Pheno` (1 = first level,
0 = second):
- `pheno_syndrome.txt`: Focal (1) vs Generalised (0)
- `pheno_DRE.txt`: DRE (1) vs no DRE (0)
- `pheno_convulsions.txt`: convulsions (1) vs none (0)
- `pheno_currentVPA.txt`: VPA exposure (1) vs none (0)

Files written by step 1 (`PRS_DIR`):
- `PRS_<gwas>_<contrast>.all_score`: per-individual PRS at 7 thresholds (`Pt_0.0001` to
  `Pt_1`). The scores do not depend on the contrast, so steps 4-6 only read the
  `_syndrome.all_score` files.
- `PRS_<gwas>_<contrast>.summary`: PRSice best-threshold result for the binary contrast.

## Configure

Copy `config.sh`, fill in absolute paths (one `KEY="value"` per line, no `${VAR}`
references, because the R and Python steps parse the file rather than source it), and pass
it as the first argument to each script. Setting `PIPELINE_CONFIG` works too.

| Key | Used by |
|---|---|
| `METADATA`, `COVARIATES` | steps 0, 1, 3-6 |
| `TARGET`, `PRSICE`, `VALID_SNPS`, `GWAS_<name>` (8) | steps 0-1 |
| `PHENO_DIR` | written by 0, read by 1 |
| `PRS_DIR` | written by 1, read by 2, 4-6 |
| `OUT_DIR` | written by 2-7 |

## Run

```
Rscript scripts/00_create_pheno_files.R config.sh
bash    scripts/01_M1_run_PRSice.sh     config.sh
Rscript scripts/02_M1_summarise_PRSice.R config.sh
Rscript scripts/03_M2_clinical_only.R   config.sh
Rscript scripts/04_M3_individual_PRS.R  config.sh
Rscript scripts/05_M4_hierarchical.R    config.sh
Rscript scripts/06_M5_stratification.R  config.sh
python3 scripts/07_FDR_correction.py    config.sh
```

If the PRSice `.all_score` files already exist, start at step 2 (or step 3).

Requirements: R with `ggplot2`, `ggpubr`, `dplyr`, `tidyr`, `reshape2`; Python 3 with
`numpy` (and `openpyxl` for the Excel copy of the FDR table); PRSice-2 for step 1.

## Outputs (in `OUT_DIR`)

| File | Step |
|---|---|
| `PRS_summary_table_noPostSurg.txt` | 2: 32 PRSice best-threshold results |
| `clinical_only_coefs.txt`, `clinical_only_forest.png` | 3 |
| `PRS_pBAG_linear_results.txt`, `PRS_pBAG_linear.pdf` | 4: 5 PRS x 7 thresholds, sorted by p |
| `hierarchical_regression_blocks.txt` | 5: cumulative and incremental R2 per block |
| `hierarchical_regression_fullmodel.txt`, `hierarchical_regression_pBAG.pdf` | 5: full model coefficients |
| `PRS_stratification_summary.txt` | 6: tertile and quartile tests per PRS |
| `PRS_interaction_summary.txt`, `PRS_stratification.pdf` | 6: PRS x clinical interaction terms |
| `FDR_correction_summary.txt`, `FDR_correction_summary.xlsx` | 7 |

Keep these out of the repository along with the inputs.
