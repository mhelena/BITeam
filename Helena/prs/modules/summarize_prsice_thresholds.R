#!/usr/bin/env Rscript
# Num_SNP and R2 at every PRSice threshold, per score, with the best Pt marked.
# Usage: Rscript summarize_prsice_thresholds.R [results_dir] [suffix]   (default: $WORKDIR/results, no suffix)
suppressPackageStartupMessages({
    library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
RES_DIR <- if (length(args) >= 1) args[1] else
           file.path(Sys.getenv("WORKDIR", "/re_gecip/<your_gecip>/<project>/<cohort>"), "results")
SUFFIX  <- if (length(args) >= 2) args[2] else ""

cat("Results dir : ", RES_DIR, "\n")
cat("File pattern: PRS_<pheno>", SUFFIX, ".prsice/.summary\n\n", sep = "")

files <- Sys.glob(file.path(RES_DIR, sprintf("PRS_*%s.prsice", SUFFIX)))
if (!length(files)) stop("No .prsice files found in ", RES_DIR)

extract_pheno <- function(f) sub(sprintf(".*PRS_(.+)%s\\.prsice$", SUFFIX), "\\1", f)

long <- do.call(rbind, lapply(files, function(f) {
    pheno <- extract_pheno(f)
    x <- read.table(f, header = TRUE, sep = "", stringsAsFactors = FALSE)
    cols_needed <- c("Threshold","Num_SNP","R2","P")
    have <- intersect(cols_needed, names(x))
    if (length(have) < 4) {
        cat("Skipping (unexpected columns):", f, "\n"); return(NULL) }
    out <- data.frame(phenotype  = pheno,
                      Threshold  = x$Threshold,
                      Num_SNP    = as.integer(x$Num_SNP),
                      R2         = as.numeric(x$R2),
                      P          = as.numeric(x$P),
                      stringsAsFactors = FALSE)
    summ_f <- sub("\\.prsice$", ".summary", f)
    if (file.exists(summ_f)) {
        s <- read.table(summ_f, header = TRUE, sep = "", stringsAsFactors = FALSE)
        out$is_best <- abs(out$Threshold - s$Threshold[1]) < 1e-9
    } else {
        out$is_best <- FALSE
    }
    out
}))

if (is.null(long) || !nrow(long)) stop("No parseable rows found.")

long <- long[order(long$phenotype, long$Threshold), ]

wide <- reshape(long[, c("phenotype","Threshold","Num_SNP")],
                idvar = "Threshold", timevar = "phenotype", direction = "wide")
names(wide) <- sub("^Num_SNP\\.", "", names(wide))
wide <- wide[order(wide$Threshold), ]

best_rows <- long[long$is_best, c("phenotype","Threshold","Num_SNP","R2","P")]
best_rows$Threshold <- format(best_rows$Threshold, scientific = FALSE)

out_wide <- file.path(RES_DIR, "snp_count_by_threshold_wide.tsv")
out_long <- file.path(RES_DIR, "snp_count_by_threshold_long.tsv")
out_best <- file.path(RES_DIR, "best_threshold_summary.tsv")

write.table(wide, out_wide, sep = "\t", row.names = FALSE, quote = FALSE)
write.table(long, out_long, sep = "\t", row.names = FALSE, quote = FALSE)
write.table(best_rows, out_best, sep = "\t", row.names = FALSE, quote = FALSE)

cat("Wrote:\n  ", out_wide, "\n  ", out_long, "\n  ", out_best, "\n\n", sep = "")

cat("Num_SNP per phenotype per threshold:\n")
print(wide, row.names = FALSE)
cat("\nBest Pt (from .summary files):\n")
print(best_rows, row.names = FALSE)

long$Threshold_lab <- format(long$Threshold, scientific = TRUE, digits = 2)
long$Threshold_lab <- factor(long$Threshold_lab,
                             levels = unique(long$Threshold_lab[order(long$Threshold)]))

# Num_SNP is drawn on the R2 axis, rescaled, with its own axis on the right
scale_factor <- with(long, max(R2, na.rm = TRUE) / max(Num_SNP, na.rm = TRUE))

g <- ggplot(long, aes(x = Threshold_lab, group = phenotype)) +
    geom_col(aes(y = R2, fill = is_best), width = 0.7) +
    geom_line(aes(y = Num_SNP * scale_factor), colour = "grey30", linewidth = 0.6) +
    geom_point(aes(y = Num_SNP * scale_factor), colour = "grey30", size = 1.4) +
    geom_text(aes(y = Num_SNP * scale_factor,
                  label = format(Num_SNP, big.mark = ",")),
              vjust = -0.6, size = 2.8, colour = "grey20") +
    scale_fill_manual(values = c(`FALSE` = "#B0BEC5", `TRUE` = "#EF6C00"),
                      labels = c("threshold", "best Pt"), name = NULL) +
    scale_y_continuous(name = expression("R"^2),
                       sec.axis = sec_axis(~ . / scale_factor,
                                            name = "N SNPs",
                                            labels = function(x) format(x, big.mark = ","))) +
    facet_wrap(~ phenotype, scales = "free_y", ncol = 2) +
    labs(x = "P-value threshold (base GWAS)",
         title = "PRSice threshold sweep: SNP count and R2 per Pt",
         subtitle = "Orange bar = best Pt from .summary. Line = N SNPs (right axis).") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          legend.position = "top",
          panel.grid.minor = element_blank(),
          strip.text = element_text(face = "bold"))

out_pdf <- file.path(RES_DIR, "snp_count_by_threshold.pdf")
out_png <- file.path(RES_DIR, "snp_count_by_threshold.png")
ggsave(out_pdf, g, width = 10, height = 7)
ggsave(out_png, g, width = 10, height = 7, dpi = 200)
cat("\nWrote:\n  ", out_pdf, "\n  ", out_png, "\n", sep = "")
