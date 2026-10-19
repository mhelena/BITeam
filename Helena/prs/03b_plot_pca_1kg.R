#!/usr/bin/env Rscript
# Plot the joint PCA against 1000 Genomes to check the samples sit in the EUR cluster.
# Usage: Rscript 03b_plot_pca_1kg.R [WORKDIR]   (in RStudio, set the defaults below)
suppressPackageStartupMessages({
    library(ggplot2)
})

args     <- commandArgs(trailingOnly = TRUE)
WORKDIR  <- if (length(args) >= 1) args[1] else Sys.getenv("WORKDIR", "/re_gecip/<your_gecip>/<project>/<cohort>")
COHORT   <- Sys.getenv("COHORT", basename(WORKDIR))
KG_PANEL <- Sys.getenv("KG_PANEL", "/public_data_resources/1000-genomes/integrated_call_samples_v3.20130502.ALL.panel")

eig   <- read.table(file.path(WORKDIR, "pca/pca_target_KG.eigenvec"), header = TRUE)
panel <- read.table(KG_PANEL, header = TRUE, sep = "\t")

eig$pop <- panel$super_pop[match(eig$IID, panel$sample)]
eig$pop[is.na(eig$pop)] <- "study"

cols <- c(AFR = "#332288", AMR = "#88CCEE", EAS = "#117733",
          EUR = "#CC6677", SAS = "#DDCC77", study = "#000000")

p <- ggplot(eig, aes(PC1, PC2, colour = pop)) +
    geom_point(data = subset(eig, pop != "study"), alpha = 0.5, size = 1.4) +
    geom_point(data = subset(eig, pop == "study"), alpha = 0.85, size = 1.4,
               shape = 18) +
    scale_colour_manual(values = cols) +
    labs(title = sprintf("%s samples projected onto 1000 Genomes", COHORT),
         x = "PC1", y = "PC2", colour = NULL) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "right")

out_png <- file.path(WORKDIR, "pca/ancestry_pca_vs_1KG.png")
out_pdf <- file.path(WORKDIR, "pca/ancestry_pca_vs_1KG.pdf")
ggsave(out_png, p, width = 7.5, height = 5.5, dpi = 300)
ggsave(out_pdf, p, width = 7.5, height = 5.5)
cat("Wrote", out_png, "\n")
