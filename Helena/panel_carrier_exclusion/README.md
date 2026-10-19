# Panel carrier exclusion

Removes likely monogenic cases from a pool of candidate controls in the GEL Research
Environment. It takes a PanelApp gene panel, finds variants in panel genes that have HIGH
or MODERATE VEP impact and are absent from gnomAD v3.1.2 (using the aggV2 VEP 109
functional annotation), and drops every pool sample carrying one of them.

## Setup

Edit `config.sh` (project folder, panel id and version, pool list) and export its path:

    export PIPELINE_CONFIG=/re_gecip/<your_gecip>/<project>/scripts/panel_carrier_exclusion/config.sh

The pool list (`POOL_FILE`) is one platekey per line.

## Run order

1. `source $PIPELINE_CONFIG && Rscript panel_get_genes.R` writes the panel gene symbols.
2. `bash panel_make_bed.sh` builds a merged CDS BED for those genes from GENCODE v38.
   Needed only when `PANEL_MODE=bed`. If `panel_ensembl.txt` (one ENSG per line) is in
   `WORK_DIR` it is used instead of the symbols.
3. `./panel_submit.sh` submits `panel_1_sites.sh` (array, one task per VEP chunk) and then
   `panel_2_carriers.sh`, which writes the carrier list and `FINAL_OUT`.
4. Optional: `source $PIPELINE_CONFIG && Rscript panel_3_final.R` when the pool was built
   from several lists (`POOL_SUBSETS`). It reports survival per list and writes a PLINK
   `--keep` file. Needs Rlabkey if any list holds participant ids.

`PANEL_MODE` decides how variants are tied to the panel. `bed` restricts the scan to
`PANEL_BED` and keeps any HIGH/MODERATE consequence there; `symbol` scans everything and
keeps rows whose VEP SYMBOL is in `panel_genes.txt`. Both apply the same impact and gnomAD
filters. Swap in your own BED (for example one that leaves out biallelic-only genes) by
pointing `PANEL_BED` at it.

## Annotation dump (optional)

`./panel_annot_submit.sh` writes every panel-gene variant (by SYMBOL, no impact or gnomAD
filter) with the pool AC, AN, AF and het/hom counts to
`${ANNOT_DIR}/_merged/panel_annot_with_pool_ac.tsv`. Useful for seeing which variants drive
the carrier counts. Stages: `1_sites` (array), `2pre_regions`, `2a_pool_ac_chunk` (array),
`2b_pool_ac_concat`, `3_merge`.

## Outputs (under `WORK_DIR`)

- `panel_sites_all.tsv`: CHROM, POS, REF, ALT, SYMBOL, Consequence, IMPACT, gnomADg_AF
- `panel_carriers.txt`: pool platekeys to exclude
- `final/controls.txt`: pool minus carriers

Jobs write `.done` or `.empty` flags per chunk, so a resubmission skips finished chunks.
Delete the flags to force a rerun.
