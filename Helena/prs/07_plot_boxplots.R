#!/usr/bin/env Rscript
# Boxplots of each PRS by cohort, and case vs control contrasts from PRS_z ~ cohort + PC1-4 + Sex.
# Usage: Rscript 07_plot_boxplots.R [WORKDIR]   (in RStudio, set the defaults below)
suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
})

args       <- commandArgs(trailingOnly = TRUE)
WORKDIR    <- if (length(args) >= 1) args[1] else Sys.getenv("WORKDIR", "/re_gecip/<your_gecip>/<project>/<cohort>")
PRS_LIST   <- strsplit(Sys.getenv("PRS_LIST", "epilepsy intelligence mvage"), "[[:space:]]+")[[1]]
CASE_LABEL <- Sys.getenv("CASE_LABEL", "Cases")
RES        <- file.path(WORKDIR, "results")

read_iid <- function(f) {
    if (!file.exists(f) || file.size(f) == 0) return(character())
    read.table(f, colClasses = "character")[, 2]
}
cases <- read_iid(file.path(WORKDIR, "merged/case_samples.txt"))
epi   <- read_iid(file.path(WORKDIR, "merged/epi_samples.txt"))
hc    <- read_iid(file.path(WORKDIR, "merged/hc_samples.txt"))
dee   <- read_iid(file.path(WORKDIR, "merged/dee_samples.txt"))

label_cohort <- function(iid) {
    ifelse(iid %in% cases, CASE_LABEL,
    ifelse(iid %in% epi,   "GEL Epi",
    ifelse(iid %in% hc,    "HC",
    ifelse(iid %in% dee,   "DEE", NA_character_))))
}

cov <- read.table(file.path(WORKDIR, "cov.txt"), header = TRUE,
                  colClasses = c(IID = "character"))

# PRS at the threshold PRSice picked (.summary), taken from .all.score
load_best_prs <- function(phe) {
    s <- read.table(file.path(RES, sprintf("PRS_%s.summary", phe)), header = TRUE)
    f <- file.path(RES, sprintf("PRS_%s.all.score", phe))
    if (!file.exists(f)) f <- file.path(RES, sprintf("PRS_%s.all_score", phe))
    x <- read.table(f, header = TRUE, check.names = FALSE, colClasses = c(IID = "character"))
    pt_cols <- setdiff(names(x), c("FID", "IID"))
    pt_vals <- as.numeric(sub("^Pt_", "", pt_cols))
    best    <- pt_cols[which.min(abs(pt_vals - s$Threshold[1]))]
    list(d = data.frame(IID = x$IID, PRS = x[[best]]),
         pt = s$Threshold[1], nsnp = s$Num_SNP[1])
}

fills <- setNames(c("#332288", "#CC6677", "#DDCC77", "#3E8E7E"),
                  c(CASE_LABEL, "GEL Epi", "HC", "DEE"))

for (phe in PRS_LIST) {
    if (!file.exists(file.path(RES, sprintf("PRS_%s.summary", phe)))) {
        message("skip ", phe, ": no PRSice output yet"); next
    }
    b <- load_best_prs(phe)
    d <- merge(b$d, cov, by = "IID")
    d$cohort <- droplevels(factor(label_cohort(d$IID),
                                  levels = c(CASE_LABEL, "GEL Epi", "HC", "DEE")))
    d <- d[!is.na(d$cohort), ]
    d$PRS_z <- as.numeric(scale(d$PRS))

    contrasts <- character()
    for (ref in setdiff(levels(d$cohort), CASE_LABEL)) {
        d$cohort_ref <- relevel(d$cohort, ref = ref)
        fit <- lm(PRS_z ~ cohort_ref + PC1 + PC2 + PC3 + PC4 + factor(Sex), data = d)
        co  <- summary(fit)$coefficients[paste0("cohort_ref", CASE_LABEL), ]
        contrasts <- c(contrasts, sprintf("%s vs %s: beta = %.3f (SE %.3f), P = %.3g",
                                          CASE_LABEL, ref, co[1], co[2], co[4]))
    }

    cat("\n", phe, ": best Pt ", b$pt, ", ", b$nsnp, " SNPs\n", sep = "")
    print(table(d$cohort))
    cat(paste0("  ", contrasts), sep = "\n")

    d$PRS_adj <- as.numeric(scale(resid(
        lm(PRS ~ PC1 + PC2 + PC3 + PC4 + factor(Sex), data = d))))

    p <- ggplot(d, aes(cohort, PRS_adj, fill = cohort)) +
        geom_boxplot(alpha = 0.8, outlier.size = 0.6) +
        scale_fill_manual(values = fills) +
        labs(title = sprintf("%s PRS (best Pt %s, %s SNPs)",
                             tools::toTitleCase(phe), b$pt, b$nsnp),
             subtitle = paste(contrasts, collapse = "\n"),
             x = NULL, y = "PRS (z-score, adjusted for PC1-4 and sex)") +
        theme_minimal(base_size = 12) +
        theme(legend.position = "none",
              plot.title = element_text(face = "bold"),
              plot.subtitle = element_text(size = 9, colour = "grey30"))

    ggsave(file.path(RES, sprintf("PRS_%s_boxplot.png", phe)), p, width = 7.5, height = 5.5, dpi = 300)
    ggsave(file.path(RES, sprintf("PRS_%s_boxplot.pdf", phe)), p, width = 7.5, height = 5.5)
    cat("  wrote results/PRS_", phe, "_boxplot.{png,pdf}\n", sep = "")
}
