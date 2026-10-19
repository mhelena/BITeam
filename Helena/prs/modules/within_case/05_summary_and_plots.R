#!/usr/bin/env Rscript
# Summary table of the within-case PRSice runs, boxplots for the binary outcome, onset scatter plots.
# Usage: Rscript 05_summary_and_plots.R   (settings from config.sh via the environment, or the defaults below)
suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
})

WORKDIR    <- Sys.getenv("WORKDIR", "/re_gecip/<your_gecip>/<project>/<cohort>")
WORK       <- Sys.getenv("WITHIN_DIR", file.path(WORKDIR, "within_case"))
OUTCOME    <- Sys.getenv("OUTCOME", "outcome")
BASES      <- strsplit(Sys.getenv("WITHIN_PRS_LIST", "bmi epilepsy intelligence longevity mvage"), "[[:space:]]+")[[1]]
RES_BIN    <- file.path(WORK, "results_binary")
RES_ONSET  <- file.path(WORK, "results_onset")
HAS_ONSET  <- file.exists(file.path(WORK, "pheno_onset.txt"))

read_prs_best <- function(res_dir, base) {
    score_f <- file.path(res_dir, sprintf("PRS_%s.all.score", base))
    summ_f  <- file.path(res_dir, sprintf("PRS_%s.summary",   base))
    if (!file.exists(score_f) || !file.exists(summ_f)) return(NULL)
    s <- read.table(summ_f, header = TRUE, sep = "", stringsAsFactors = FALSE)
    best_pt <- s$Threshold[1]
    r2      <- if ("PRS.R2"  %in% names(s)) s$PRS.R2[1]  else NA
    nsnp    <- if ("Num_SNP" %in% names(s)) s$Num_SNP[1] else NA
    P_run   <- if ("P"       %in% names(s)) s$P[1]       else NA
    coef    <- if ("Coefficient" %in% names(s)) s$Coefficient[1] else NA
    se      <- if ("Standard.Error" %in% names(s)) s$Standard.Error[1] else NA
    x <- read.table(score_f, header = TRUE, check.names = FALSE,
                    stringsAsFactors = FALSE)
    if (!"IID" %in% names(x)) names(x)[1:2] <- c("FID","IID")
    pt_cols <- setdiff(names(x), c("FID","IID"))
    pt_nums <- suppressWarnings(as.numeric(sub("^Pt_", "", pt_cols)))
    idx <- which.min(abs(pt_nums - best_pt))
    best_col <- pt_cols[idx]
    list(scores = data.frame(IID = x$IID, PRS = x[[best_col]],
                             stringsAsFactors = FALSE),
         best_pt = best_pt, r2 = r2, nsnp = nsnp, P = P_run,
         beta = coef, se = se)
}

fmt_p <- function(x) {
    if (is.na(x)) return("NA")
    if (x < 1e-4) return(format(x, scientific = TRUE, digits = 2))
    formatC(x, format = "fg", digits = 2)
}

cov <- read.table(file.path(WORK, "cov.txt"), header = TRUE)
oc  <- read.table(file.path(WORK, "pheno_binary.txt"), header = TRUE)
oc$Y_lab <- factor(ifelse(oc$Outcome == 2, "Yes", "No"), levels = c("No","Yes"))

bin_rows <- list()
for (base in BASES) {
    info <- read_prs_best(RES_BIN, base)
    if (is.null(info)) { message("MISS ", OUTCOME, " ", base); next }
    d <- info$scores %>%
        inner_join(oc[, c("IID","Y_lab")], by = "IID") %>%
        inner_join(cov, by = "IID") %>%
        filter(!is.na(PRS), !is.na(Y_lab))
    d$PRS_z <- as.numeric(scale(d$PRS))
    fit <- lm(PRS_z ~ Y_lab + PC1 + PC2 + PC3 + PC4 + factor(Sex), data = d)
    co <- summary(fit)$coefficients
    term <- grep("^Y_lab", rownames(co), value = TRUE)[1]
    beta_z <- co[term, "Estimate"]; se_z <- co[term,"Std. Error"]
    P_z <- co[term, "Pr(>|t|)"]
    n_yes <- sum(d$Y_lab == "Yes"); n_no <- sum(d$Y_lab == "No")
    bin_rows[[length(bin_rows)+1]] <- data.frame(
        Phenotype = OUTCOME, PRS = base,
        best_Pt = info$best_pt, N_SNP = info$nsnp,
        PRS_R2 = signif(info$r2, 3),
        N_case = n_yes, N_ctrl = n_no,
        beta_z = signif(beta_z, 3), SE_z = signif(se_z, 3),
        P_adj  = signif(P_z, 3),
        stringsAsFactors = FALSE)

    fit_adj <- lm(PRS ~ PC1 + PC2 + PC3 + PC4 + factor(Sex), data = d)
    d$PRS_adj <- as.numeric(scale(residuals(fit_adj)))
    means_df <- aggregate(PRS_adj ~ Y_lab, data = d, FUN = mean)
    means_df$label <- sprintf("%.2f", means_df$PRS_adj)
    y_max <- max(d$PRS_adj); by <- y_max + 0.5
    labs <- setNames(sprintf("%s-%s (n=%d)", OUTCOME, levels(d$Y_lab), c(n_no, n_yes)),
                     levels(d$Y_lab))
    g <- ggplot(d, aes(Y_lab, PRS_adj, color = Y_lab)) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
        geom_jitter(width = 0.22, alpha = 0.55, size = 0.9) +
        geom_boxplot(fill = NA, outlier.shape = NA, width = 0.5, linewidth = 0.6) +
        stat_summary(fun = mean, geom = "point", shape = 23, size = 3.5,
                     fill = "white", stroke = 1, color = "black",
                     show.legend = FALSE) +
        geom_text(data = means_df, inherit.aes = FALSE,
                  aes(x = Y_lab, y = PRS_adj, label = label),
                  vjust = -1.3, hjust = -0.15, size = 3.2, fontface = "bold") +
        annotate("segment", x = 1, xend = 2, y = by, yend = by, color = "grey25") +
        annotate("segment", x = 1, xend = 1, y = by, yend = by - 0.12, color = "grey25") +
        annotate("segment", x = 2, xend = 2, y = by, yend = by - 0.12, color = "grey25") +
        annotate("text", x = 1.5, y = by + 0.18,
                 label = fmt_p(P_z), size = 3.6, color = "grey15") +
        scale_color_manual(values = c("No" = "#0072B2", "Yes" = "#E69F00"),
                           labels = labs, name = "Phenotype:") +
        labs(title = sprintf("PRS - %s (within cases, %s)   best Pt=%s, R2=%.3g, SNPs=%s",
                              tools::toTitleCase(base), OUTCOME,
                              format(info$best_pt), info$r2, format(info$nsnp)),
             x = NULL,
             y = sprintf("PRS - %s (adjusted z)", tools::toTitleCase(base))) +
        theme_minimal(base_size = 12) +
        theme(legend.position = "right", panel.grid.minor = element_blank(),
              plot.title = element_text(face = "bold", size = 11))
    ggsave(file.path(WORK, sprintf("PRS_%s_%s_boxplot.png", base, OUTCOME)),
           g, width = 7, height = 5, dpi = 300)
    ggsave(file.path(WORK, sprintf("PRS_%s_%s_boxplot.pdf", base, OUTCOME)),
           g, width = 7, height = 5)
    cat("Wrote", OUTCOME, "boxplot:", base, "\n")
}
bin_tab <- do.call(rbind, bin_rows)

combined <- data.frame(Phenotype = bin_tab$Phenotype, PRS = bin_tab$PRS,
                       best_Pt = bin_tab$best_Pt, N_SNP = bin_tab$N_SNP,
                       PRS_R2 = bin_tab$PRS_R2,
                       N_case_or_total = bin_tab$N_case,
                       N_ctrl_or_NA    = bin_tab$N_ctrl,
                       beta = bin_tab$beta_z, SE = bin_tab$SE_z,
                       P = bin_tab$P_adj, stringsAsFactors = FALSE)

if (HAS_ONSET) {
    onset <- read.table(file.path(WORK, "pheno_onset.txt"), header = TRUE)

    # The primary fit uses every case; the onset cutoffs are sensitivity analyses
    # for a long right tail.
    q3   <- quantile(onset$onset_months, 0.75, na.rm = TRUE)
    iqr  <- IQR(onset$onset_months, na.rm = TRUE)
    iqr_thr <- unname(q3 + 1.5 * iqr)

    cutoffs <- list(
        full = list(label = "No exclusion",       keep = function(x) rep(TRUE, length(x)), thr = NA),
        gt18 = list(label = "Onset <= 18 mo",     keep = function(x) x <= 18,              thr = 18),
        gt12 = list(label = "Onset <= 12 mo",     keep = function(x) x <= 12,              thr = 12),
        iqr  = list(label = sprintf("Onset <= %.1f mo (Q3+1.5*IQR)", iqr_thr),
                    keep = function(x) x <= iqr_thr, thr = iqr_thr)
    )

    fit_one <- function(d, keep_fn) {
        d2 <- d[keep_fn(d$onset_months), ]
        if (nrow(d2) < 30) return(list(n = nrow(d2), effect = NA, se = NA, p = NA))
        d2$PRS_z <- as.numeric(scale(d2$PRS))
        fit <- lm(PRS_z ~ onset_months + PC1 + PC2 + PC3 + PC4 + factor(Sex),
                  data = d2)
        s <- summary(fit)$coefficients
        list(n = nrow(d2),
             effect = s["onset_months","Estimate"],
             se     = s["onset_months","Std. Error"],
             p      = s["onset_months","Pr(>|t|)"])
    }

    onset_rows <- list()
    for (base in BASES) {
        info <- read_prs_best(RES_ONSET, base)
        if (is.null(info)) { message("MISS ONSET ", base); next }
        d <- info$scores %>%
            inner_join(onset[, c("IID","onset_months")], by = "IID") %>%
            inner_join(cov, by = "IID") %>%
            filter(!is.na(PRS), !is.na(onset_months))
        if (nrow(d) < 30) next

        rows <- lapply(names(cutoffs), function(k) {
            r <- fit_one(d, cutoffs[[k]]$keep)
            data.frame(PRS = base, cutoff = k, cutoff_label = cutoffs[[k]]$label,
                       n_kept = r$n, n_excluded = nrow(d) - r$n,
                       effect_per_month = r$effect, SE = r$se, P = r$p,
                       stringsAsFactors = FALSE)
        })
        onset_rows[[base]] <- do.call(rbind, rows)

        fit_adj <- lm(PRS ~ PC1 + PC2 + PC3 + PC4 + factor(Sex), data = d)
        d$PRS_adj <- as.numeric(scale(residuals(fit_adj)))
        d$excluded_12 <- d$onset_months > 12

        ann <- onset_rows[[base]]
        primary_txt <- sprintf("Primary (all cases, n=%d): beta=%+.3f, p=%s",
                               ann$n_kept[1], ann$effect_per_month[1],
                               fmt_p(ann$P[1]))
        sens_txt <- paste(
            sprintf("<=18 mo (n=%d): beta=%+.3f, p=%s",
                    ann$n_kept[2], ann$effect_per_month[2], fmt_p(ann$P[2])),
            sprintf("<=12 mo (n=%d): beta=%+.3f, p=%s",
                    ann$n_kept[3], ann$effect_per_month[3], fmt_p(ann$P[3])),
            sprintf("<=IQR   (n=%d): beta=%+.3f, p=%s",
                    ann$n_kept[4], ann$effect_per_month[4], fmt_p(ann$P[4])),
            sep = "\n")

        g <- ggplot(d, aes(onset_months, PRS_adj)) +
            geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
            geom_vline(xintercept = 12, linetype = "dotted", color = "grey60") +
            geom_vline(xintercept = 18, linetype = "dotted", color = "grey60") +
            geom_point(color = "#0072B2", alpha = 0.7, size = 1.6) +
            geom_smooth(method = "lm", se = TRUE,
                        color = "black", linewidth = 0.9) +
            geom_smooth(data = subset(d, !excluded_12),
                        method = "lm", se = FALSE,
                        color = "grey55", linewidth = 0.6, linetype = "dashed") +
            annotate("label", x = -Inf, y = Inf, label = primary_txt,
                     hjust = -0.02, vjust = 1.1, size = 3.4, fontface = "bold",
                     fill = "white", label.size = 0.3) +
            annotate("label", x = Inf, y = Inf, label = sens_txt,
                     hjust = 1.02, vjust = 1.05, size = 2.9, family = "mono",
                     fill = "white", label.size = 0.3) +
            labs(title = sprintf("PRS - %s vs age at seizure onset (cases)",
                                  tools::toTitleCase(base)),
                 subtitle = sprintf("Primary fit: all cases (solid black). Reference: <=12 mo only (grey dashed). best Pt=%s, R2=%.3g, SNPs=%s",
                                     format(info$best_pt), info$r2, format(info$nsnp)),
                 x = "Age at seizure onset (months)",
                 y = sprintf("PRS - %s (adjusted z)", tools::toTitleCase(base))) +
            theme_minimal(base_size = 12) +
            theme(plot.title    = element_text(face = "bold"),
                  plot.subtitle = element_text(size = 9, color = "grey30"),
                  panel.grid.minor = element_blank(),
                  legend.position = "none")

        ggsave(file.path(WORK, sprintf("PRS_%s_ONSET_scatter.png", base)),
               g, width = 8, height = 6, dpi = 300)
        ggsave(file.path(WORK, sprintf("PRS_%s_ONSET_scatter.pdf", base)),
               g, width = 8, height = 6)
        cat("Wrote onset scatter:", base, "\n")
    }

    onset_sens_tab <- do.call(rbind, onset_rows)
    write.table(onset_sens_tab,
                file.path(WORK, "PRS_by_onset_sensitivity_table.tsv"),
                sep = "\t", row.names = FALSE, quote = FALSE)
    cat("\nOnset sensitivity table:\n")
    print(onset_sens_tab, row.names = FALSE)

    onset_summary <- onset_sens_tab %>%
        filter(cutoff == "full") %>%
        transmute(Phenotype = "Onset (months)", PRS,
                  best_Pt = NA, N_SNP = NA, PRS_R2 = NA,
                  N = n_kept,
                  beta_per_month = signif(effect_per_month, 3),
                  SE = signif(SE, 3),
                  P_adj = signif(P, 3))
    for (i in seq_len(nrow(onset_summary))) {
        info <- read_prs_best(RES_ONSET, onset_summary$PRS[i])
        if (!is.null(info)) {
            onset_summary$best_Pt[i] <- info$best_pt
            onset_summary$N_SNP[i]   <- info$nsnp
            onset_summary$PRS_R2[i]  <- signif(info$r2, 3)
        }
    }

    combined <- rbind(combined,
        data.frame(Phenotype = onset_summary$Phenotype, PRS = onset_summary$PRS,
                   best_Pt = onset_summary$best_Pt, N_SNP = onset_summary$N_SNP,
                   PRS_R2 = onset_summary$PRS_R2,
                   N_case_or_total = onset_summary$N,
                   N_ctrl_or_NA    = NA,
                   beta = onset_summary$beta_per_month, SE = onset_summary$SE,
                   P = onset_summary$P_adj, stringsAsFactors = FALSE))
}

write.table(combined, file.path(WORK, "within_case_summary_table.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("\nWithin-case summary:\n")
print(combined, row.names = FALSE)
cat("\nWrote", file.path(WORK, "within_case_summary_table.tsv"), "\n")
