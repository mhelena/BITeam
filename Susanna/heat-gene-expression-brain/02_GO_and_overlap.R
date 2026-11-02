
library(dplyr)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(UpSetR)

datasets <- c("passive_heat_1h", "heatstroke", "neurons_1h", "neurons_3h", "neurons_24h", "organoids")
de <- lapply(datasets, function(n) read.delim(paste0("DE/DE_", n, ".txt")))
names(de) <- datasets

up   <- lapply(de, function(x) x$gene[x$adj.P.Val < 0.05 & x$logFC > 0])
down <- lapply(de, function(x) x$gene[x$adj.P.Val < 0.05 & x$logFC < 0])
sig  <- lapply(de, function(x) x$gene[x$adj.P.Val < 0.05])

############################################################
# 1. GENES SHARED BY THE FOUR SYSTEMS
############################################################

four <- c("passive_heat_1h", "heatstroke", "neurons_1h", "organoids")
shared <- Reduce(intersect, sig[four])
length(shared)
length(Reduce(intersect, down[four])) 
length(Reduce(intersect, up[four]))

# same thing if neurons = all 3 timepoints together
neur_down <- unique(unlist(down[c("neurons_1h","neurons_3h","neurons_24h")]))
neur_up   <- unique(unlist(up[c("neurons_1h","neurons_3h","neurons_24h")]))
length(Reduce(intersect, list(down$passive_heat_1h, down$heatstroke, neur_down, down$organoids)))
length(Reduce(intersect, list(up$passive_heat_1h, up$heatstroke, neur_up, up$organoids)))

# pairwise overlap and direction agreement (Figure 1B)
for (i in 1:3) for (j in (i+1):4) {
  a <- de[[four[i]]]; b <- de[[four[j]]]
  m <- inner_join(a, b, by = "gene") %>% filter(adj.P.Val.x < 0.05, adj.P.Val.y < 0.05)
  cat(four[i], four[j], "shared:", nrow(m),
      " same direction:", round(100 * mean(sign(m$logFC.x) == sign(m$logFC.y))), "%\n")
}

pdf("SupplFig10_upset_4systems.pdf", width = 7, height = 4)
upset(fromList(down[four]), nsets = 4, nintersects = 15, order.by = "freq",
      mainbar.y.label = "Genes in intersection", sets.x.label = "Set size")
upset(fromList(up[four]), nsets = 4, nintersects = 15, order.by = "freq",
      mainbar.y.label = "Genes in intersection", sets.x.label = "Set size")
dev.off()

############################################################
# 2. GO BP
############################################################

run_go <- function(genes) {
  eg <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  ego <- enrichGO(eg$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP",
                  pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05,
                  readable = TRUE)
  filter(as.data.frame(ego), p.adjust < 0.05)
}

go_up   <- lapply(up, run_go)
go_down <- lapply(down, run_go)

sapply(go_up, nrow)
sapply(go_down, nrow)

saveRDS(list(up = go_up, down = go_down), "GO_per_dataset.rds")

ggplot(head(go_down$heatstroke, 8),
       aes(Count / as.numeric(sub(".*/", "", GeneRatio)), reorder(Description, -p.adjust))) +
  geom_point(aes(size = Count, color = -log10(p.adjust))) +
  scale_color_viridis_c(option = "plasma", direction = -1) +
  theme_bw() +
  labs(title = "Clinical heatstroke - downregulated",
       x = "Gene ratio (genes in term / genes tested)", y = "GO BP term")

pool_neurons <- function(go) {
  bind_rows(go$neurons_1h, go$neurons_3h, go$neurons_24h) %>%
    group_by(Description) %>% slice_min(p.adjust, n = 1, with_ties = FALSE) %>%
    ungroup() %>% arrange(p.adjust)
}
neur_down_go <- pool_neurons(go_down)
neur_up_go   <- pool_neurons(go_up) 
neur_terms <- neur_down_go$Description

shared_terms <- Reduce(intersect, list(go_down$passive_heat_1h$Description,
                                       go_down$heatstroke$Description,
                                       neur_terms,
                                       go_down$organoids$Description))
length(shared_terms)

# same with neurons at 1h only - gives the same 72
length(Reduce(intersect, list(go_down$passive_heat_1h$Description, go_down$heatstroke$Description,
                              go_down$neurons_1h$Description, go_down$organoids$Description)))
length(Reduce(intersect, list(go_up$passive_heat_1h$Description, go_up$heatstroke$Description,
                              neur_up_go$Description, go_up$organoids$Description)))

go_rank <- function(go, terms) {
  go <- go[order(go$p.adjust, go$Description, method = "radix"), ]
  match(terms, go$Description)
}
rank_tab <- data.frame(Description = shared_terms,
                       passive_heat = go_rank(go_down$passive_heat_1h, shared_terms),
                       heatstroke   = go_rank(go_down$heatstroke, shared_terms),
                       neurons      = go_rank(neur_down_go, shared_terms),
                       organoids    = go_rank(go_down$organoids, shared_terms))
head(rank_tab[order(apply(rank_tab[, -1], 1, min)), ])
write.table(rank_tab, "shared_72_terms_ranks.txt", row.names = F, sep = '\t', quote = F)

# heat shock / folding terms in the upregulated lists (Supplementary Table 2)
hs_terms <- c("response to topologically incorrect protein",
              "response to unfolded protein",
              "'de novo' post-translational protein folding",
              "protein refolding",
              "chaperone-mediated protein folding",
              "response to temperature stimulus",
              "'de novo' protein folding",
              "response to heat",
              "chaperone cofactor-dependent protein refolding",
              "protein folding",
              "cellular response to heat",
              "regulation of cellular response to heat",
              "protein stabilization",
              "cellular response to topologically incorrect protein",
              "cellular response to unfolded protein",
              "chaperone-mediated protein complex assembly")
data.frame(term         = hs_terms,
           passive_heat = go_rank(go_up$passive_heat_1h, hs_terms),
           neurons      = go_rank(neur_up_go, hs_terms),
           organoids    = go_rank(go_up$organoids, hs_terms),
           heatstroke   = go_rank(go_up$heatstroke, hs_terms))

fig2 <- bind_rows(
  mutate(filter(go_down$passive_heat_1h, Description %in% shared_terms), system = "Passive heat PBMCs"),
  mutate(filter(neur_down_go, Description %in% shared_terms), system = "Primary cortical neurons"),
  mutate(filter(go_down$organoids, Description %in% shared_terms), system = "Cortical organoids"),
  mutate(filter(go_down$heatstroke, Description %in% shared_terms), system = "Clinical heatstroke PBMCs"))
fig2$ratio <- fig2$Count / as.numeric(sub(".*/", "", fig2$GeneRatio))
top25 <- fig2 %>% group_by(Description) %>% summarise(m = median(-log10(p.adjust))) %>%
  arrange(desc(m)) %>% head(25)
fig2 <- filter(fig2, Description %in% top25$Description)
fig2$Description <- factor(fig2$Description, levels = rev(top25$Description))
fig2$system <- factor(fig2$system, levels = c("Passive heat PBMCs", "Primary cortical neurons","Cortical organoids", "Clinical heatstroke PBMCs"))

ggplot(fig2, aes(system, Description)) +
  geom_point(aes(size = ratio, color = -log10(p.adjust))) +
  scale_color_gradient(low = "#d6e6f5", high = "darkred", limits = c(-log10(0.05), NA)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  labs(x = "Heat-stress system", y = "GO Biological Process term",size = "Gene ratio", color = "-log10 (BH-adjusted P)")
ggsave("Fig2.pdf", width = 170, height = 130, units = "mm")

