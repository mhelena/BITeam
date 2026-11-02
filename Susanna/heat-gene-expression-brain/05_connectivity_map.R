
library(dplyr)
library(signatureSearch)
library(ExperimentHub)
library(rhdf5)
library(org.Hs.eg.db)
library(httr)

# LINCS reference from ExperimentHub
eh <- ExperimentHub()
q <- query(eh, c("signatureSearchData", "lincs"))
lincs <- eh[[names(q)[q$title == "lincs"][1]]]
lincs_genes <- as.integer(h5read(lincs, "rownames"))

to_entrez <- function(g) {
  e <- unname(na.omit(mapIds(org.Hs.eg.db, unique(g), "ENTREZID", "SYMBOL")))
  as.character(intersect(e, lincs_genes))
}

run_lincs <- function(up, down) {
  qs <- qSig(query = list(upset = to_entrez(up), downset = to_entrez(down)),
             gess_method = "LINCS", refdb = lincs)
  as.data.frame(result(gess_lincs(qs, sortby = "NCS", tau = TRUE, workers = 4)))
}

maxq <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA)
  qh <- quantile(x, 0.67); ql <- quantile(x, 0.33)
  unname(if (abs(qh) >= abs(ql)) qh else ql)
}
summarise_tau <- function(tb) {
  tb %>% group_by(pert) %>%
    summarise(n_cell = n(), tau = maxq(Tau), best_tau = suppressWarnings(min(Tau, na.rm = TRUE)),
              fdr = median(WTCS_FDR, na.rm = TRUE), moa = dplyr::first(MOAss)) %>%
    arrange(tau)
}

############################################################
# 1. HEATSTROKE SIGNATURE
############################################################

hs <- read.delim("DE/DE_heatstroke.txt")
sig_up   <- head(hs$gene[hs$adj.P.Val < 0.05 & hs$logFC >  0.5], 150)
sig_down <- head(hs$gene[hs$adj.P.Val < 0.05 & hs$logFC < -0.5], 150)

hs_lincs <- run_lincs(sig_up, sig_down)
hs_tau <- summarise_tau(hs_lincs)

reversers <- filter(hs_tau, tau <= -90, fdr < 0.05, n_cell >= 2)
write.table(reversers, "SuppTable6_heatstroke_reversers.txt", row.names = F, sep = '\t', quote = F)

############################################################
# 2. TSC / mTORC1 activation signature
############################################################

enrichr_set <- function(library, term) {
  txt <- content(GET(paste0("https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=", library)), "text", encoding = "UTF-8")
  lines <- grep(term, strsplit(txt, "\n")[[1]], value = TRUE)
  g <- unlist(lapply(lines, function(l) strsplit(l, "\t")[[1]][-(1:2)]))
  unique(toupper(g[g != ""]))
}
tsc_up   <- enrichr_set("MSigDB_Hallmark_2020", "^mTORC1 Signaling\t")
tsc_down <- enrichr_set("KEGG_2021_Human", "^(Autophagy|Lysosome)\t")
both <- intersect(tsc_up, tsc_down)
tsc_up <- setdiff(tsc_up, both); tsc_down <- setdiff(tsc_down, both)

tsc_tau <- summarise_tau(run_lincs(tsc_up, tsc_down))
filter(tsc_tau, tolower(pert) %in% c("sirolimus", "temsirolimus"))

############################################################
# 3. GRANULE vs PURKINJE SIGNATURE 
############################################################

div <- read.delim("divergent_genes.txt")
gp <- div %>% group_by(gene) %>%
  summarise(n_datasets = n_distinct(dataset), logFC = median(logFC)) %>%
  filter(n_datasets >= 2) %>%
  arrange(desc(n_datasets), gene)

gp_up   <- head(gp$gene[gp$logFC > 0], 150)
gp_down <- head(gp$gene[gp$logFC < 0], 150)

gp_tau <- summarise_tau(run_lincs(gp_up, gp_down))
gp_reversers <- filter(gp_tau, tau <= -90, fdr < 0.05, n_cell >= 2)
write.table(gp_reversers, "SuppTable7_granule_purkinje_reversers.txt", row.names = F, sep = '\t', quote = F)
