# PRS within a case-only cohort: workshop

This workshop computes polygenic risk scores (PRS) in a cohort of people with epilepsy **without controls**, and compares the PRS between **groups defined by genotype**, for example people whose epilepsy is explained by a single rare variant against people with no such variant. Changing which group plays the "case" and which the "comparison" lets you ask several questions with the same scores.

Everything runs inside RStudio. plink2 and PRSice are called from the notebooks, so you never need a separate terminal.

## Setup (once)

1. Open `prs_case_only.Rproj` in RStudio. This sets the working folder, so all paths work.
2. Install the R packages if you do not have them:
   ```r
   install.packages(c("dplyr", "tidyr", "readr", "ggplot2", "stringr", "patchwork", "forcats",
                      "purrr", "scales", "rlang", "shiny", "rmarkdown", "knitr"))
   ```
3. Install **plink2** (https://www.cog-genomics.org/plink/2.0/) and **PRSice-2** (https://choishingwan.github.io/PRSice/).
4. Open `config.R` and set `PLINK2` and `PRSICE` to the full path of each program (only needed if they are not on your PATH).
5. Open `00_make_toy_data.R` and click **Source**. It builds the practice cohort in `data/toy/` in about a minute.

## The workshop

Open each notebook in turn and run it **one chunk at a time** (the green Run arrow, or Ctrl/Cmd + Shift + Enter), reading the output before moving on. **Your turn** chunks are exercises, with a suggested answer hidden in a comment.

| Step | File | What you do |
|---|---|---|
| 0 | `00_make_toy_data.R` | Build the simulated cohort (Source it once) |
| 1 | `01_explore_cohort.Rmd` | Tidy the clinical data, summarise genetic findings per person, define the genotype groups |
| 2 | `02_genotype_qc.Rmd` | Missingness, allele frequency, Hardy-Weinberg, heterozygosity, relatedness |
| 3 | `03_ancestry_pcs.Rmd` | Principal components, and why ancestry matters when groups are defined by genotype |
| 4 | `04_prs_scores.Rmd` | Clean the GWAS, compute the scores with PRSice (scores only, no threshold search) |
| 5 | `05_within_cohort_contrasts.Rmd` | Pre-specified comparisons, negative control, PC adjustment, threshold sensitivity, power |
| 6 | `06_app.R` | Click **Run App**: pick any two groups and see the comparison update live |

Results, plots and logs go to `results/`. `R/functions.R` holds the helper functions (`run_tool`, `compare_groups`, `plot_groups`, `forest_plot`, `score_col`).

## Using your own cohort

Edit the data block in `config.R`:

- `GENO`: PLINK binary genotypes (`.bed/.bim/.fam`), autosomes, one ancestry group or with PCs.
- `CLINICAL`: one row per person. Notebook 1 renames the columns, so adapt the `rename()` line to your column names.
- `FINDINGS`: one row per reported rare variant with columns `IID, gene, variant, consequence, acmg`.
- `GWAS`: base summary statistics with columns `SNP, CHR, BP, A1 (effect allele), A2, FRQ, BETA, SE, P`. Rename yours to match in notebook 4, and use a GWAS that did **not** include your cohort.
- Gene lists and `PT_MAIN` (the main P value threshold): set them **before** looking at any PRS.

## Key ideas

- **No controls, so no threshold search.** PRSice normally picks the P value threshold that best separates cases from controls. We ask only for the scores (`--no-regress`), fix the main threshold in advance and show the others as a sensitivity check.
- **Liability threshold.** Everyone crossed the same threshold to get epilepsy, so people with a strong rare variant needed less polygenic risk. Inside a case-only cohort, a monogenic group should have a lower PRS.
- **Write the hypotheses and their directions first.** Every extra comparison is another test.
- **Negative-control PRS.** A trait unrelated to epilepsy (height here) should show no difference. If it does, look for ancestry or batch effects.
- **Adjust for PCs.** Genotype groups can differ in ancestry (founder variants), and ancestry shifts the PRS.

---

## Notes for the workshop leader

The toy cohort is simulated with a known truth (`data/toy/INSTRUCTOR_truth.tsv` and `INSTRUCTOR_qc_problems.tsv`; you may want to delete these before sharing the folder):

- **1,201 genotyped samples**: 1,200 people plus one duplicate. The clinical table has 1,210 rows: the 1,200 people plus 10 with clinical data but no genotypes. The duplicate has genotypes but no clinical record.
- **QC should remove 13 samples**: 6 with high missingness, 3 contaminated (excess heterozygosity), 3 relatives and the duplicate (`_R`). The 3 contaminated samples look "related" to each other if heterozygosity is not checked first, a point made in notebook 2.
- **Simulated truth**: high-penetrance carriers have the lowest polygenic liability, then incomplete-penetrance genes (DEPDC5, NPRL2, NPRL3), then VUS carriers and people with no variant (these two are the same: a VUS has no effect in the simulation). Missense variants in high-penetrance genes are slightly milder than PTVs.
- **NPRL3** is a founder variant in the smaller ancestry cluster, which has a slightly lower average PRS. H3 is underestimated without PCs (about 0.32 SD) and recovered with them (about 0.41 SD).
- **Expected results** at P < 0.1 (seed 2026):

  | Comparison | Epilepsy PRS | Height PRS |
  |---|---|---|
  | H1 no variant vs monogenic | +0.46 SD, P about 1e-10 | about 0 |
  | H2 VUS only vs monogenic | +0.39 SD, P about 0.02 | about 0 |
  | H3 incomplete vs high penetrance | +0.41 SD, P about 0.002 | about 0 |
  | H4 missense vs PTV | +0.41 SD, P about 0.04 | +0.34 SD, P 0.07 (chance, small groups) |

  The H4 height result is a useful discussion point: with groups of 40 to 70, a sizeable difference can appear by chance, and the interval shows it.
- The effect is stronger at stricter thresholds in this toy (the simulated GWAS is small, so loose thresholds add mostly noise). Real epilepsy PRS effects are smaller, typically 0.1 to 0.3 SD.
- Tested end to end on macOS with R 4.4.1, plink2 v2.0.0-a.7 and PRSice 2.3.5.
