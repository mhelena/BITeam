#!/usr/bin/env Rscript
# Write the HGNC symbols of a PanelApp panel at a pinned version to PANEL_GENES.
# Usage: source "$PIPELINE_CONFIG" && Rscript panel_get_genes.R
suppressPackageStartupMessages({
  library(httr); library(jsonlite)
})

env <- function(k) {
  v <- Sys.getenv(k)
  if (!nzchar(v)) stop(k, " is not set; source config.sh first")
  v
}
PANEL_ID <- env("PANEL_ID")
VER      <- env("PANEL_VERSION")
CONF     <- env("PANEL_CONFIDENCE")
out      <- env("PANEL_GENES")
dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)

url <- sprintf("https://panelapp.genomicsengland.co.uk/api/v1/panels/%s/?version=%s",
               PANEL_ID, VER)
cat(sprintf("[%s] fetching %s\n", format(Sys.time()), url))

r <- GET(url, timeout(60))
stop_for_status(r)
panel <- fromJSON(content(r, "text", encoding = "UTF-8"), simplifyVector = FALSE)

genes_raw <- panel$genes
cat(sprintf("  panel '%s' v%s: %d entities total\n",
            panel$name, panel$version, length(genes_raw)))

keep_levels <- if (CONF == "Green") "3" else if (CONF == "Amber") c("2","3") else c("1","2","3")
sel <- vapply(genes_raw, function(g) as.character(g$confidence_level) %in% keep_levels,
              logical(1))
genes <- vapply(genes_raw[sel], function(g) g$gene_data$gene_symbol, character(1))
genes <- sort(unique(genes))

writeLines(genes, out)
cat(sprintf("  kept %d genes at confidence >= %s -> %s\n",
            length(genes), CONF, out))
