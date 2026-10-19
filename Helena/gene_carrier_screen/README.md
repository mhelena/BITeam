# Gene carrier screen (GEL RE)

Finds GEL rare-disease probands from disease categories unrelated to a phenotype who carry
missense variants in one gene that are absent from gnomAD, then builds a manual review
sheet, an HPO table and an anonymised variant table for Airlock export. Written for SCN1A
in epilepsy (the control definition from PMID 37006128) but the gene, region, transcript
and phenotype terms are all set in `config.sh`.

## Inputs

- aggV2 genomic_data and VEP 109 functional_annotation (paths in `config.sh`).
- LabKey access through Rlabkey and `~/.netrc`.
- `WHITELIST`: disease categories unrelated to the phenotype, one per line. Carriers are
  kept only if their GEL disease (specific disease, subgroup or group) or an HPO term
  matches a line. For epilepsy this is Supplementary Table 3 of PMID 37006128. The file is
  not in the repository.
- `EXCLUDE_IDS` (optional): platekeys or participant IDs to drop, one per line.

## Run order

Copy `config.sh`, fill in the paths, then:

```bash
export PIPELINE_CONFIG=/absolute/path/to/config.sh
source "$PIPELINE_CONFIG"
module load R/4.2.0

Rscript 01_build_pool.R          # hg38 probands without the phenotype, any ancestry
bash submit.sh scan              # LSF: 02_sites.sh then 03_carriers.sh
Rscript 04_filter_controls.R     # whitelist and EUR >= 0.8, four output lists
Rscript 05_hpo_per_sample.R      # HPO long/wide tables, age filter
bash submit.sh annotate          # LSF: 06_annotate_variants.sh
Rscript 07_review_sheet.R        # one row per carrier x variant
Rscript 08_airlock_table.R       # anonymised table in AIRLOCK_DIR
```

Steps 4 and 5 accept another ID list as their first argument.

## Filters

- Sites: `SYMBOL == GENE`, consequence contains `missense_variant` on any transcript, and
  `gnomADg_AF` empty. The gnomAD release is the one in the GEL annotation (v3.1.2 for
  VEP 109), so a few sites may be present in newer gnomAD.
- Annotation: one row per variant on the canonical protein-coding transcript (or
  `TRANSCRIPT` if set), excluding `NMD_transcript_variant`. Carriers whose variant has no
  such row are flagged `<no_canonical_missense>` in the review sheet.
- Controls: the whitelist is applied with an intersect, so only matching participants are
  kept.

## Outputs (in OUT_DIR)

| File | Contents |
|---|---|
| `pool.txt` | candidate pool, platekeys |
| `carriers.txt`, `carrier_variants.tsv` | carriers and their variants |
| `carriers_{allanc,whitelist,eur08,whitelist_eur08}.txt` | filtered lists |
| `*_hpo_long.tsv`, `*_hpo_wide.tsv` | HPO terms per sample |
| `variants_annotated.tsv` | canonical annotation per variant |
| `review_sheet_<GENE>_<date>.tsv` | manual review sheet |
| `for_airlock/<GENE>_controls_airlock_variants.tsv` | anonymised table for export |
| `<GENE>_controls_sample_id_mapping.tsv` | Sample ID to platekey map, never export |

Everything except the anonymised table contains participant identifiers and stays in the RE.
