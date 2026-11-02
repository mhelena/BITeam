
library(data.table)
library(dplyr)
library(ggplot2)


gtex_file <- "GTEx_Analysis_2025-08-22_v11_RNASeQCv2.4.3_gene_median_tpm.gct.gz"

gtex <- fread(gtex_file, skip = 2, data.table = FALSE)

brain_cols <- grep("Brain", colnames(gtex), value = TRUE)
gtex <- gtex[, c("Description", brain_cols)]

gtex$mean <- rowMeans(gtex[, brain_cols])
gtex <- gtex[order(-gtex$mean), ]
gtex <- gtex[!duplicated(gtex$Description), ]

brain <- as.matrix(gtex[, brain_cols])
rownames(brain) <- gtex$Description
colnames(brain) <- gsub("Brain_|_basal_ganglia|_cervical_c-1", "", colnames(brain))

brain_log <- log2(brain + 1)
brain_log <- brain_log[rowMeans(brain_log) > 0.5, ] 
nrow(brain_log)

############################################################
# 1. AUC WITHIN EACH REGION
# the cerebellum has higher expression overall, so we do not compare absolute values: we rank all genes inside each region and ask where the heat genes fall (AUC = Mann-Whitney U / (n_deg * n_background)) 0.5 = no enrichment
############################################################

ranks <- apply(brain_log, 2, rank)
pct   <- ranks / nrow(ranks)

region_auc <- function(genes) {
  g  <- intersect(genes, rownames(ranks))
  n1 <- length(g)
  n2 <- nrow(ranks) - n1
  (colSums(ranks[g, , drop = FALSE]) - n1 * (n1 + 1) / 2) / (n1 * n2)
}

datasets <- c("passive_heat_1h", "neurons_1h", "organoids", "neurons_3h","heatstroke", "neurons_24h")
de <- lapply(datasets, function(n) read.delim(paste0("DE/DE_", n, ".txt")))
names(de) <- datasets

direction <- "down"
#direction <- "up"
genes <- lapply(de, function(x) {
  if (direction == "up") x$gene[x$adj.P.Val < 0.05 & x$logFC > 0]
  else x$gene[x$adj.P.Val < 0.05 & x$logFC < 0]
})

auc <- sapply(genes, region_auc)
round(auc, 3)
apply(-auc, 2, rank)[c("Cerebellum", "Cerebellar_Hemisphere"), ]

plot_df <- data.frame(region = rep(rownames(auc), ncol(auc)),
                      dataset = rep(colnames(auc), each = nrow(auc)),
                      auc = as.vector(auc))
region_order <- names(sort(rowMeans(auc)))
plot_df$region <- factor(plot_df$region, levels = region_order)
plot_df$dataset <- factor(plot_df$dataset, levels = datasets)

ggplot(plot_df, aes(auc, region)) +
  geom_vline(xintercept = 0.5, linetype = "dashed", colour = "grey45") +
  geom_point(colour = "#5B9BD5") +
  facet_wrap(~ dataset, ncol = 2) +
  theme_minimal() +
  labs(x = "Within-region AUC", y = "Brain region")
ggsave(paste0("Fig_regions_", direction, ".pdf"), width = 170, height = 200, units = "mm")

############################################################
# 2. PAIRWISE REGION COMPARISONS 
############################################################

pairwise <- list()
for (n in datasets) {
  g <- intersect(genes[[n]], rownames(pct))
  pr <- combn(colnames(pct), 2)
  p <- apply(pr, 2, function(r) wilcox.test(pct[g, r[1]], pct[g, r[2]], paired = TRUE)$p.value)
  pairwise[[n]] <- data.frame(dataset = n, region1 = pr[1, ], region2 = pr[2, ],
                              diff = auc[pr[1, ], n] - auc[pr[2, ], n],
                              p_adj = p.adjust(p, "BH"))
}
pairwise <- bind_rows(pairwise)
pairwise$stars <- as.character(cut(pairwise$p_adj, c(0, 0.001, 0.01, 0.05, 1), labels = c("***","**","*",""), include.lowest = TRUE))

pairwise %>%
  filter(region1 == "Cerebellum" | region2 == "Cerebellum", p_adj < 0.05) %>%
  mutate(cb_lower = ifelse(region1 == "Cerebellum", diff < 0, diff > 0)) %>%
  group_by(dataset) %>%
  summarise(below = sum(cb_lower), above = sum(!cb_lower))

write.table(pairwise, paste0("region_pairwise_", direction, ".txt"), row.names = F, sep = '\t', quote = F)

d <- filter(pairwise, dataset == "heatstroke")
d <- bind_rows(d, mutate(d, region1 = region2, region2 = d$region1, diff = -diff))
d <- filter(d, match(region2, region_order) > match(region1, region_order))
d$region1 <- factor(d$region1, levels = region_order)
d$region2 <- factor(d$region2, levels = region_order)
ggplot(d, aes(region2, region1)) +
  geom_tile(aes(fill = ifelse(p_adj < 0.05, -diff, NA)), colour = "white") +
  geom_text(aes(label = stars), size = 2.5) +
  scale_fill_gradient2(low = "#2471a3", mid = "white", high = "#c0392b", na.value = "grey92",name = "Difference in AUC\n(column minus row)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Brain region", y = "Brain region", title = "Clinical heatstroke")
