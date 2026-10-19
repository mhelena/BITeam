# LDSC genetic correlation

Pairwise genetic correlations (rg) between GWAS summary statistics with LDSC,
plus a heatmap and forest plot. Every pair of traits in the trait table is run.

## Needs

- LDSC (github.com/bulik/ldsc) and its python environment
- European HapMap3 LD scores and weights (`eur_w_ld_chr` style, hg19)
- R with ggplot2 for the plots

## Trait table

One row per GWAS, tab separated. Two examples are in `traits/`: `traits_ds.tsv`
(the Dravet PRS traits) and `traits_brainage.tsv` (pBAG project).

| column | meaning |
|---|---|
| name | short name used for files |
| label | name shown on plots |
| file | path relative to `GWAS_DIR`, can be gzipped, tab or whitespace separated |
| snp_col | SNP column |
| snp_type | `rsid`, or `chrbp` for chr:pos IDs (GRCh37), mapped to rsIDs through the HapMap3 LD score files |
| a1_col / a2_col | effect / other allele |
| effect_col | effect column |
| effect_type | `beta`, or `or` (converted to log OR) |
| p_col | p value column |
| n | a column name, a fixed N, or cases/controls such as `29944/52538` |

## Running

```
cp config.sh my_config.sh     # edit paths
PIPELINE_CONFIG=$PWD/my_config.sh bash run_all.sh
```

Steps, if you want to run them one at a time:

1. `01_prep_sumstats.py` writes `prep/<name>.txt` (SNP A1 A2 BETA P [N]) and `prep/traits.list`
2. `02_munge.sh` runs munge_sumstats.py into `munged/`
3. `03_rg_all_pairs.sh` runs `ldsc.py --rg` for every pair into `rg/`
4. `04_collect.sh` builds `rg_summary.tsv` (trait1 trait2 rg se z p h2_obs1 h2_obs2)
5. `05_plot_rg.R` writes `rg_heatmap.png`, `rg_forest.png` and `rg_plots.pdf`

`PLOT_FOCUS` limits the forest plot to pairs involving the listed traits
(e.g. `brainage_gm,brainage_wm,brainage_gwm`). The heatmap always shows all pairs.

## Notes

- rg is noisy when a GWAS has h2 z < 4 (see the LDSC FAQ). Check
  the h2 columns before reading too much into a pair.
- N for case-control GWAS: give cases/controls or an effective N column, not both.
- No summary statistics or results go in the repo.
