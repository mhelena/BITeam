#!/usr/bin/env Rscript
# Heatmap and forest plot of the pairwise rg table.
# usage: Rscript 05_plot_rg.R rg_summary.tsv traits.tsv out_dir [focus_traits] [rg_limit]

suppressPackageStartupMessages(library(ggplot2))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("usage: 05_plot_rg.R rg_summary.tsv traits.tsv out_dir [focus] [limit]")
rg_file <- args[1]
traits  <- read.delim(args[2], stringsAsFactors = FALSE, comment.char = "#")
out_dir <- args[3]
focus   <- if (length(args) >= 4 && nzchar(args[4])) strsplit(args[4], ",")[[1]] else character(0)
lim     <- if (length(args) >= 5 && nzchar(args[5])) as.numeric(args[5]) else 1

rg <- read.delim(rg_file, stringsAsFactors = FALSE)
for (v in c("rg", "se", "z", "p")) rg[[v]] <- suppressWarnings(as.numeric(rg[[v]]))
rg <- rg[!is.na(rg$rg), ]

lab <- setNames(traits$label, traits$name)
lvl <- traits$label

stars <- function(p) {
  ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", ifelse(p < 0.1, ".", ""))))
}
rg$sig <- stars(rg$p)

# both triangles so the heatmap reads the same from either axis
both <- rbind(
  data.frame(a = lab[rg$trait1], b = lab[rg$trait2], rg = rg$rg, sig = rg$sig),
  data.frame(a = lab[rg$trait2], b = lab[rg$trait1], rg = rg$rg, sig = rg$sig)
)
both$a <- factor(both$a, levels = lvl)
both$b <- factor(both$b, levels = rev(lvl))
# values outside the limit get squished to the end colour rather than greyed out
both$fill <- pmax(pmin(both$rg, lim), -lim)

p1 <- ggplot(both, aes(a, b, fill = fill)) +
  geom_tile(colour = "white", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.2f%s", rg, sig)), size = 3.8) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-lim, lim), name = expression(r[g])) +
  labs(x = NULL, y = NULL, title = "LDSC genetic correlations",
       subtitle = "*** p<0.001, ** p<0.01, * p<0.05, . p<0.1") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        panel.grid = element_blank(),
        plot.title = element_text(face = "bold"))

fr <- rg
if (length(focus)) fr <- fr[fr$trait1 %in% focus | fr$trait2 %in% focus, ]
fr$pair <- paste(lab[fr$trait1], "vs", lab[fr$trait2])
fr$pair <- factor(fr$pair, levels = fr$pair[order(fr$rg)])
fr$lo <- fr$rg - 1.96 * fr$se
fr$hi <- fr$rg + 1.96 * fr$se
fr$grp <- factor(ifelse(fr$p < 0.05, "p < 0.05", ifelse(fr$p < 0.1, "p < 0.1", "n.s.")),
                 levels = c("p < 0.05", "p < 0.1", "n.s."))

p2 <- ggplot(fr, aes(rg, pair, colour = grp)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.25, linewidth = 0.7, orientation = "y") +
  geom_point(size = 3) +
  scale_colour_manual(values = c("p < 0.05" = "#B2182B", "p < 0.1" = "#F4A582", "n.s." = "#636363"),
                      name = NULL, drop = FALSE) +
  labs(x = expression(r[g] ~ "(95% CI)"), y = NULL, title = "LDSC genetic correlations") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

n <- length(lvl)
ggsave(file.path(out_dir, "rg_heatmap.png"), p1, width = 2 + 1.1 * n, height = 1.5 + 0.9 * n, dpi = 300)
ggsave(file.path(out_dir, "rg_forest.png"), p2, width = 8, height = 2 + 0.3 * nrow(fr), dpi = 300)

pdf(file.path(out_dir, "rg_plots.pdf"), width = 2 + 1.1 * n, height = max(6, 2 + 0.3 * nrow(fr)))
print(p1)
print(p2)
invisible(dev.off())
