
library(GEOquery)
library(limma)
library(dplyr)
library(tibble)
library(data.table)
library(AnnotationDbi)
library(hgu133plus2.db)
library(org.Hs.eg.db)
library(ggplot2)
library(ggrepel)
library(patchwork)


############################################################
# 1. PASSIVE HEAT EXPOSURE (SAUNA) PBMC - GSE90763
############################################################

gset <- getGEO("GSE90763", GSEMatrix = TRUE)[[1]]
expr <- exprs(gset)   # already log2 normalized (RMA)
meta <- pData(gset)

meta$condition <- case_when(
  grepl("before", meta$source_name_ch1) ~ "baseline",
  grepl("15 min", meta$source_name_ch1) ~ "heat_15min",
  grepl("1 hr", meta$source_name_ch1) ~ "heat_1h")
meta$condition <- factor(meta$condition, levels = c("baseline","heat_15min","heat_1h"))
meta$subject <- factor(meta$`subject identifier:ch1`)

# sex and age of the volunteers (one row per subject)
subj <- meta[!duplicated(meta$subject), ]
table(subj$`gender:ch1`)
range(as.numeric(subj$`age:ch1`))

symbols <- mapIds(hgu133plus2.db, keys = rownames(expr), column = "SYMBOL",
                  keytype = "PROBEID", multiVals = "first")
expr_df <- as.data.frame(expr)
expr_df$Gene <- symbols
expr_df <- filter(expr_df, !is.na(Gene) & Gene != "")

expr_df$var <- apply(expr_df[, colnames(expr)], 1, var)
expr_gene <- expr_df %>%
  group_by(Gene) %>%
  slice_max(order_by = var, n = 1, with_ties = FALSE) %>%
  ungroup()
expr_mat <- as.matrix(expr_gene[, colnames(expr)])
rownames(expr_mat) <- expr_gene$Gene

design <- model.matrix(~ subject + condition, data = meta)
fit <- eBayes(lmFit(expr_mat, design))

res_15 <- topTable(fit, coef = "conditionheat_15min", number = Inf)
res_1h <- topTable(fit, coef = "conditionheat_1h", number = Inf)

sum(res_15$adj.P.Val < 0.05)
sum(res_1h$adj.P.Val < 0.05)

passive_15 <- rownames_to_column(res_15, "gene")
passive_1h <- rownames_to_column(res_1h, "gene")

############################################################
# 2. CLINICAL HEATSTROKE PBMC - GSE317880
############################################################

gset <- getGEO("GSE317880", GSEMatrix = TRUE)[[1]]
expr <- exprs(gset)
meta <- pData(gset)

# patients on arrival (before cooling) vs heat-exposed controls
keep <- meta$source_name_ch1 %in% c("Control", "Heat Stroke T0")
meta <- meta[keep, ]
expr <- expr[, keep]

meta$group <- factor(ifelse(meta$source_name_ch1 == "Control", "Control", "HS"),
                     levels = c("Control","HS"))
meta$age <- as.numeric(meta$`age:ch1`)
meta$sex <- factor(meta$`gender:ch1`)
table(meta$group, meta$sex)
tapply(meta$age, meta$group, range)

# adjust for age and sex
design <- model.matrix(~ group + age + sex, data = meta)
fit <- eBayes(lmFit(expr, design))
tt <- topTable(fit, coef = "groupHS", number = Inf, sort.by = "none")

# Clariom D probes are RefSeq ids (e.g. NM_000016.6_at) -> gene symbol
refseq <- sub("\\.[0-9]+$", "", sub("_at$", "", rownames(tt)))
tt$gene <- mapIds(org.Hs.eg.db, keys = refseq, column = "SYMBOL",
                  keytype = "REFSEQ", multiVals = "first")

heatstroke <- tt %>%
  filter(!is.na(gene) & gene != "") %>%
  group_by(gene) %>%
  slice_min(P.Value, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(P.Value)

sum(heatstroke$adj.P.Val < 0.05)

############################################################
# 3. PRIMARY NEURONS - GSE132447
############################################################

cuff <- fread("neurons/GSE132447_gene_exp.diff.gz", data.table = FALSE)
cuff$sample_1 <- tolower(cuff$sample_1)
cuff$sample_2 <- tolower(cuff$sample_2)
cuff$logFC <- as.numeric(cuff$`log2(fold_change)`)
cuff$adj.P.Val <- as.numeric(cuff$q_value)

neurons <- list()
for (tp in c("hs1hr", "hs3hr", "hs24hr")) {

  x <- filter(cuff, sample_1 == "untreated", sample_2 == tp, status == "OK",
              gene != "-", is.finite(logFC))
  x$refseq <- sub(",.*", "", x$gene)
  x$gene <- mapIds(org.Hs.eg.db, keys = x$refseq, column = "SYMBOL",
                   keytype = "REFSEQ", multiVals = "first")

  neurons[[tp]] <- x %>%
    filter(!is.na(gene)) %>%
    group_by(gene) %>%
    slice_min(adj.P.Val, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    dplyr::select(gene, logFC, adj.P.Val)
}
sapply(neurons, function(x) sum(x$adj.P.Val < 0.05))

############################################################
# 4. CORTICAL ORGANOIDS - GSE149797
############################################################

tpm <- fread("neurons/GSE149797_all.tpm.txt.gz", data.table = FALSE)

ctrl_cols <- c("control_1", "control_2", "control_3")
hs_cols   <- c("HS_control_1", "HS_control_2", "HS_control_3")

tpm <- tpm[tpm$Symbol != "" & !duplicated(tpm$Symbol), ]
tpm_mat <- as.matrix(tpm[, c(ctrl_cols, hs_cols)])
rownames(tpm_mat) <- tpm$Symbol

tpm_log <- log2(tpm_mat + 1)

keep <- rowMeans(tpm_log[, ctrl_cols]) > 1 | rowMeans(tpm_log[, hs_cols]) > 1
tpm_log <- tpm_log[keep, ]

group <- factor(c(rep("ctrl", 3), rep("hs", 3)))
design <- model.matrix(~ group)
fit <- eBayes(lmFit(tpm_log, design))
organoids <- rownames_to_column(topTable(fit, coef = "grouphs", number = Inf), "gene")

sum(organoids$adj.P.Val < 0.05)

############################################################
# 5. FINAL
############################################################

de <- list(
  passive_heat_15min = passive_15,
  passive_heat_1h    = passive_1h,
  heatstroke         = heatstroke,
  neurons_1h         = neurons$hs1hr,
  neurons_3h         = neurons$hs3hr,
  neurons_24h        = neurons$hs24hr,
  organoids          = organoids)

dir.create("DE", showWarnings = FALSE)
for (n in names(de)) {
  write.table(de[[n]][, c("gene","logFC","adj.P.Val")], paste0("DE/DE_", n, ".txt"),
              row.names = F, sep = '\t', quote = F)
}

deg_counts <- data.frame(
  dataset = names(de),
  up   = sapply(de, function(x) sum(x$adj.P.Val < 0.05 & x$logFC > 0)),
  down = sapply(de, function(x) sum(x$adj.P.Val < 0.05 & x$logFC < 0)))
deg_counts$total <- deg_counts$up + deg_counts$down
deg_counts

plots <- list()
for (n in c("passive_heat_1h", "heatstroke", "organoids",
            "neurons_1h", "neurons_3h", "neurons_24h")) {

  res <- de[[n]]
  res$group <- "Not sig"
  res$group[res$adj.P.Val < 0.05 & res$logFC > 0] <- "Up"
  res$group[res$adj.P.Val < 0.05 & res$logFC < 0] <- "Down"
  top_up <- res %>% filter(group == "Up") %>% arrange(adj.P.Val, desc(abs(logFC))) %>% head(5)
  top_dn <- res %>% filter(group == "Down") %>% arrange(adj.P.Val, desc(abs(logFC))) %>% head(5)
  res$label <- ifelse(res$gene %in% c(top_up$gene, top_dn$gene), res$gene, NA)
  res$logFC <- pmax(pmin(res$logFC, 8), -8)

  plots[[n]] <- ggplot(res, aes(logFC, -log10(adj.P.Val))) +
    geom_point(aes(color = group), size = 0.6) +
    scale_color_manual(values = c("Up" = "#c0392b", "Down" = "#2471a3", "Not sig" = "lightgrey")) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_text_repel(aes(label = label), size = 2.5, max.overlaps = Inf) +
    xlim(-8, 8) +
    theme_bw() +
    theme(legend.position = "none") +
    labs(title = n, x = "log2 fold change", y = "-log10 (FDR-adjusted P)")
}

p <- wrap_plots(plots, ncol = 2, byrow = FALSE)
ggsave("SupplFig_volcano.pdf", p, width = 170, height = 185, units = "mm")
