
library(data.table)
library(dplyr)
library(ggplot2)
library(lme4)
library(emmeans)
library(limma)
library(edgeR)
library(rhdf5)
library(Matrix)
library(GEOquery)
library(clusterProfiler)
library(org.Hs.eg.db)

gse97930_file  <- "GSE97930_CerebellarHem_snDrop-seq_UMI_Count_Matrix_08-01-2017.txt.gz"
gse246183_file <- "GSE246183_human_ctype.h5ad"

datasets <- c("passive_heat_1h", "neurons_1h", "organoids", "neurons_3h","heatstroke", "neurons_24h")
de <- lapply(datasets, function(n) read.delim(paste0("DE/DE_", n, ".txt")))
names(de) <- datasets
up   <- lapply(de, function(x) x$gene[x$adj.P.Val < 0.05 & x$logFC > 0])
down <- lapply(de, function(x) x$gene[x$adj.P.Val < 0.05 & x$logFC < 0])

############################################################
# 1. GSE97930 PSEUDOBULK
############################################################

hdr <- readLines(gzfile(gse97930_file), n = 1)
cells <- gsub('"', '', strsplit(hdr, "\t")[[1]])
raw <- fread(gse97930_file, skip = 1, header = FALSE, data.table = FALSE)
counts <- as.matrix(raw[, -1])
rownames(counts) <- gsub('"', '', raw[[1]])
colnames(counts) <- cells
rm(raw)

cluster  <- sub("_cbm[0-9]+_.*", "", cells)
celltype <- cluster
celltype[grepl("^Purk", celltype)] <- "Purk"
library_id <- regmatches(cells, regexpr("cbm[0-9]+", cells))
table(celltype)

donor_of <- c(cbm1 = "0006-YO", cbm2 = "0006-YO", cbm3 = "0007-OX", cbm4 = "0007-OX", cbm5 = "0008-CR",
              cbm6 = "0008-CR", cbm7 = "0008-CR", cbm8 = "0009-IU", cbm9 = "0009-IU")

keep_ct <- c(Gran = "Granule neurons", Purk = "Purkinje neurons", Ast_Cer = "Cerebellar astrocytes",
             Ast = "Astrocytes", Oli = "Oligodendrocytes", OPC = "Oligodendrocyte precursor cells",
             OPC_Cer = "Cerebellar OPCs", Mic = "Microglia", End = "Endothelial cells", Per = "Pericytes")

combos <- unique(data.frame(lib = library_id, ct = celltype))
combos <- combos[combos$ct %in% names(keep_ct), ]
lcpm <- function(cols) {
  x <- rowSums(counts[, cols, drop = FALSE])
  if (sum(x) < 100) return(NULL)
  log2(x / sum(x) * 1e6 + 1)
}
pb <- lapply(1:nrow(combos), function(i) lcpm(library_id == combos$lib[i] & celltype == combos$ct[i]))
combos <- combos[!sapply(pb, is.null), ]
pb <- do.call(cbind, pb)
colnames(pb) <- paste(combos$ct, combos$lib, sep = "__")
pb_all <- pb
pb <- pb[rowMeans(pb) > 0.5, ]

############################################################
# 2. AUC WITHIN EACH CELL TYPE
############################################################

ct_auc <- function(genes) {
  is_deg <- rownames(pb) %in% genes
  nd <- sum(is_deg); nb <- nrow(pb) - nd
  apply(pb, 2, function(col) { r <- rank(col); (sum(r[is_deg]) - nd * (nd + 1) / 2) / (nd * nb) })
}

direction <- "down"
#direction <- "up"
genes <- if (direction == "up") up else down

auc_res <- list()
pairwise <- list()
for (n in datasets) {
  d <- data.frame(ct = keep_ct[combos$ct], lib = combos$lib, donor = donor_of[combos$lib],
                  auc = ct_auc(genes[[n]]))
  full <- lmer(auc ~ ct + (1 | donor) + (1 | donor:ct), data = d, REML = FALSE)
  null <- lmer(auc ~ 1  + (1 | donor) + (1 | donor:ct), data = d, REML = FALSE)
  print(anova(null, full)) 
  pairwise[[n]] <- mutate(as.data.frame(pairs(emmeans(full, ~ ct), adjust = "BH")), dataset = n)
  auc_res[[n]] <- as.data.frame(emmeans(full, ~ ct)) %>%
    dplyr::select(ct, mean = emmean, se = SE) %>%
    mutate(dataset = n) %>% arrange(desc(mean))
}
auc_res <- bind_rows(auc_res)
filter(auc_res, ct %in% c("Granule neurons", "Purkinje neurons"))
write.table(bind_rows(pairwise), paste0("celltype_pairwise_", direction, ".txt"), row.names = F, sep = '\t', quote = F)

ggplot(auc_res, aes(mean, reorder(ct, mean))) +
  geom_vline(xintercept = 0.5, linetype = "dashed") +
  geom_errorbarh(aes(xmin = mean - se, xmax = mean + se), height = 0) +
  geom_point(colour = "#5B9BD5") +
  facet_wrap(~ factor(dataset, levels = datasets), ncol = 2) +
  theme_minimal() +
  labs(x = "Within-cell-type AUC", y = "Cerebellar cell type")

############################################################
# 3. TAU SPECIFICITY
############################################################
ct_mean <- sapply(names(keep_ct), function(k) rowMeans(pb_all[, combos$ct == k, drop = FALSE]))
tau <- apply(ct_mean, 1, function(x) if (max(x) <= 0) NA else sum(1 - x / max(x)) / (length(x) - 1))
identity_genes <- names(tau)[!is.na(tau) & tau >= 0.8]

############################################################
# 4. GRANULE vs PURKINJE IN GSE97930
############################################################

libs <- intersect(combos$lib[combos$ct == "Gran"], combos$lib[combos$ct == "Purk"])
gran <- sapply(libs, function(l) lcpm(library_id == l & cluster == "Gran"))
purk <- sapply(libs, function(l) {
  p <- lapply(c("Purk1", "Purk2"), function(k) lcpm(library_id == l & cluster == k))
  rowMeans(do.call(cbind, p))
})
gp <- cbind(gran, purk)
rownames(gp) <- rownames(counts)
gp <- gp[rowMeans(gp) > 0.5, ]
cell <- factor(rep(c("Gran", "Purk"), each = length(libs)), levels = c("Gran", "Purk"))
blk  <- unname(donor_of[rep(libs, 2)])
design <- model.matrix(~ cell)
dc <- duplicateCorrelation(gp, design, block = blk)
dc$consensus.correlation
fit <- eBayes(lmFit(gp, design, block = blk, correlation = dc$consensus.correlation))
gvp <- topTable(fit, coef = "cellPurk", number = Inf)
gvp$gene <- rownames(gvp)

mean(gvp$logFC > 0)

heat <- list(
  PBMC       = de$passive_heat_1h$gene[de$passive_heat_1h$adj.P.Val < 0.05],
  Heatstroke = de$heatstroke$gene[de$heatstroke$adj.P.Val < 0.05],
  Neurons    = unique(c(de$neurons_1h$gene[de$neurons_1h$adj.P.Val < 0.05],
                        de$neurons_3h$gene[de$neurons_3h$adj.P.Val < 0.05],
                        de$neurons_24h$gene[de$neurons_24h$adj.P.Val < 0.05])),
  Organoids  = de$organoids$gene[de$organoids$adj.P.Val < 0.05])

go_bp <- function(genes, universe) {
  if (length(genes) < 10) return(NULL)
  eg <- bitr(genes, "SYMBOL", "ENTREZID", org.Hs.eg.db)$ENTREZID
  un <- bitr(universe, "SYMBOL", "ENTREZID", org.Hs.eg.db)$ENTREZID
  ego <- enrichGO(eg, org.Hs.eg.db, universe = un, ont = "BP", pAdjustMethod = "BH",
                  pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE)
  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) return(NULL)
  ego <- simplify(ego, cutoff = 0.7, by = "p.adjust", select_fun = min)
  filter(as.data.frame(ego), p.adjust < 0.05)
}

divergent_go <- function(ref, fc = 0) {
  out <- list()
  for (n in names(heat)) {
    g <- intersect(heat[[n]], ref$gene)
    r <- ref[match(g, ref$gene), ]
    purk <- setdiff(r$gene[r$adj.P.Val < 0.05 & r$logFC >=  fc & r$logFC > 0], identity_genes)
    gran <- setdiff(r$gene[r$adj.P.Val < 0.05 & r$logFC <= -fc & r$logFC < 0], identity_genes)
    cat(n, " Purkinje-higher:", length(purk), " granule-higher:", length(gran), "\n")
    out[[paste(n, "Purkinje")]] <- go_bp(purk, g)
    out[[paste(n, "granule")]]  <- go_bp(gran, g)
  }
  bind_rows(out, .id = "set")
}

go_97930 <- divergent_go(gvp)
table(go_97930$set)

div <- bind_rows(lapply(names(heat), function(n) {
  r <- gvp[gvp$gene %in% heat[[n]] & gvp$adj.P.Val < 0.05, ]
  data.frame(gene = r$gene, dataset = n, logFC = r$logFC, identity = r$gene %in% identity_genes)
}))
table(div$dataset, ifelse(div$logFC > 0, "Purkinje", "granule"), div$identity)
write.table(filter(div, !identity), "divergent_genes.txt", row.names = F, sep = '\t', quote = F)

go_nofilter <- divergent_go(mutate(gvp, adj.P.Val = 0))
table(go_nofilter$set)

syn <- "synap|neurotransmit|membrane potential|ion transport|ion transmembrane|celljunction|axon|dendrit|excitab|calcium|neuron projection"
go_97930 %>% filter(grepl(syn, Description)) %>% count(set)

############################################################
# 5. REPLICATION - GSE246183
############################################################

obs_cat <- function(col) {
  codes <- h5read(gse246183_file, paste0("obs/", col))
  cats  <- h5read(gse246183_file, paste0("obs/__categories/", col))
  as.vector(cats[codes + 1])
}
ct2      <- obs_cat("sub4")       # PC = Purkinje, GC = granule
genotype <- obs_cat("genotype")
donor    <- obs_cat("batch")
donor[donor == "4092"]  <- "409"
donor[donor == "15162"] <- "1516"

# control donors with at least 20 Purkinje nuclei (6 donors)
table(donor[genotype == "Ctrl" & ct2 == "PC"])
donors <- names(which(table(donor[genotype == "Ctrl" & ct2 == "PC"]) >= 20))
use <- genotype == "Ctrl" & ct2 %in% c("PC", "GC") & donor %in% donors
group <- ifelse(use, paste(donor, ct2, sep = "_"), NA)

genes246 <- as.vector(h5read(gse246183_file, "raw/var/_index"))
indptr   <- as.vector(h5read(gse246183_file, "raw/X/indptr"))
pb246 <- matrix(0, length(genes246), length(unique(na.omit(group))),
                dimnames = list(genes246, sort(unique(na.omit(group)))))
chunk <- 20000
for (s in seq(1, length(group), by = chunk)) {
  e <- min(s + chunk - 1, length(group))
  if (all(is.na(group[s:e]))) next
  a <- indptr[s]; b <- indptr[e + 1]
  x <- h5read(gse246183_file, "raw/X/data",    start = a + 1, count = b - a)
  j <- h5read(gse246183_file, "raw/X/indices", start = a + 1, count = b - a)
  m <- sparseMatrix(i = rep(1:(e - s + 1), diff(indptr[s:(e + 1)])), j = j + 1, x = as.numeric(x),
                    dims = c(e - s + 1, length(genes246)))
  g <- group[s:e]
  for (k in unique(na.omit(g))) pb246[, k] <- pb246[, k] + colSums(m[which(g == k), , drop = FALSE])
}
colSums(pb246)

# edgeR + limma-voom, donor paired
pb246 <- rowsum(pb246, rownames(pb246))
y <- DGEList(pb246)
cell <- factor(ifelse(grepl("_PC$", colnames(pb246)), "Purkinje", "Granule"), levels = c("Granule", "Purkinje"))
dnr  <- factor(sub("_.*", "", colnames(pb246)))
y <- y[filterByExpr(y, group = cell), , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)
design <- model.matrix(~ dnr + cell)
ref246 <- topTable(eBayes(lmFit(voom(y, design), design)), coef = "cellPurkinje", number = Inf)
ref246$gene <- rownames(ref246)
nrow(ref246)
mean(ref246$logFC > 0)

go_246 <- divergent_go(ref246)
table(go_246$set)

# unique synaptic terms on each side (all 4 datasets together)
go_246 %>% filter(grepl(syn, Description)) %>% mutate(side = sub(".* ", "", set)) %>%
  distinct(side, Description) %>% count(side)

