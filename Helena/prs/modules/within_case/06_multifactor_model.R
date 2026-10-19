#!/usr/bin/env Rscript
# Joint logistic model of the binary outcome on all JOINT_PRS_LIST scores (each at its own best Pt) + PC1-4 + Sex.
# Usage: Rscript 06_multifactor_model.R   (settings from config.sh via the environment, or the defaults below)
suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
})

WORKDIR <- Sys.getenv("WORKDIR", "/re_gecip/<your_gecip>/<project>/<cohort>")
WORK    <- Sys.getenv("WITHIN_DIR", file.path(WORKDIR, "within_case"))
OUTCOME <- Sys.getenv("OUTCOME", "outcome")
BASES   <- strsplit(Sys.getenv("JOINT_PRS_LIST", "epilepsy intelligence mvage bmi"), "[[:space:]]+")[[1]]
RES_BIN <- file.path(WORK, "results_binary")
PFX     <- sprintf("within_case_%s_multifactor", OUTCOME)
NB      <- length(BASES)

read_prs_best <- function(base) {
    score_f <- file.path(RES_BIN, sprintf("PRS_%s.all.score", base))
    summ_f  <- file.path(RES_BIN, sprintf("PRS_%s.summary",   base))
    if (!file.exists(score_f) || !file.exists(summ_f))
        stop(sprintf("MISSING %s PRSice output for %s", OUTCOME, base))
    s <- read.table(summ_f, header = TRUE, sep = "", stringsAsFactors = FALSE)
    best_pt <- s$Threshold[1]
    r2      <- if ("PRS.R2"  %in% names(s)) s$PRS.R2[1]  else NA
    nsnp    <- if ("Num_SNP" %in% names(s)) s$Num_SNP[1] else NA
    x <- read.table(score_f, header = TRUE, check.names = FALSE,
                    stringsAsFactors = FALSE)
    if (!"IID" %in% names(x)) names(x)[1:2] <- c("FID","IID")
    pt_cols <- setdiff(names(x), c("FID","IID"))
    pt_nums <- suppressWarnings(as.numeric(sub("^Pt_", "", pt_cols)))
    idx <- which.min(abs(pt_nums - best_pt))
    best_col <- pt_cols[idx]
    list(scores = data.frame(IID = x$IID, PRS = x[[best_col]],
                             stringsAsFactors = FALSE),
         best_pt = best_pt, r2 = r2, nsnp = nsnp)
}

cov <- read.table(file.path(WORK, "cov.txt"), header = TRUE,
                  stringsAsFactors = FALSE)
oc <- read.table(file.path(WORK, "pheno_binary.txt"), header = TRUE,
                 stringsAsFactors = FALSE)
oc$Y01 <- ifelse(oc$Outcome == 2, 1L, ifelse(oc$Outcome == 1, 0L, NA_integer_))
oc <- oc[!is.na(oc$Y01), c("IID","Y01")]
cat(sprintf("%s samples: %d (yes=%d, no=%d)\n",
            OUTCOME, nrow(oc), sum(oc$Y01 == 1), sum(oc$Y01 == 0)))

d <- oc
best_pt_by_base <- c(); r2_by_base <- c(); nsnp_by_base <- c()
for (base in BASES) {
    info <- read_prs_best(base)
    col  <- sprintf("PRS_%s", base)
    tmp  <- info$scores
    names(tmp)[names(tmp) == "PRS"] <- col
    d <- merge(d, tmp, by = "IID", all.x = TRUE)
    best_pt_by_base[base] <- info$best_pt
    r2_by_base[base]      <- info$r2
    nsnp_by_base[base]    <- info$nsnp
}
d <- merge(d, cov, by = "IID", all.x = TRUE)

must_have <- c(sprintf("PRS_%s", BASES), "PC1","PC2","PC3","PC4","Sex","Y01")
d <- d[complete.cases(d[, must_have]), ]
cat(sprintf("After complete-case filter: n=%d (yes=%d, no=%d)\n",
            nrow(d), sum(d$Y01 == 1), sum(d$Y01 == 0)))

# z-scored within the cases so the coefficients are all per SD
for (base in BASES) {
    col <- sprintf("PRS_%s",     base)
    z   <- sprintf("PRS_%s_z",   base)
    d[[z]] <- as.numeric(scale(d[[col]]))
}

# Sex alone can carry the baseline ROC, so print the 2x2 as a check
cat("\nSex x ", OUTCOME, " (n=", nrow(d), "):\n", sep = "")
sex_tab <- table(Sex = d$Sex, Outcome = ifelse(d$Y01 == 1, "yes", "no"))
print(sex_tab)
sex_ft  <- tryCatch(fisher.test(sex_tab), error = function(e) NULL)
if (!is.null(sex_ft)) {
    cat(sprintf("Fisher's exact: OR = %.3f  (95%% CI %.3f, %.3f)   P = %.3g\n",
                as.numeric(sex_ft$estimate),
                sex_ft$conf.int[1], sex_ft$conf.int[2], sex_ft$p.value))
}

form_joint    <- as.formula(paste0(
    "Y01 ~ ",
    paste(sprintf("PRS_%s_z", BASES), collapse = " + "),
    " + PC1 + PC2 + PC3 + PC4 + factor(Sex)"))
form_base     <- as.formula("Y01 ~ PC1 + PC2 + PC3 + PC4 + factor(Sex)")
form_base_pcs <- as.formula("Y01 ~ PC1 + PC2 + PC3 + PC4")
form_base_sex <- as.formula("Y01 ~ factor(Sex)")

fit_joint    <- glm(form_joint,    family = binomial, data = d)
fit_base     <- glm(form_base,     family = binomial, data = d)
fit_base_pcs <- glm(form_base_pcs, family = binomial, data = d)
fit_base_sex <- glm(form_base_sex, family = binomial, data = d)

# LR test of the PRS block over PCs + Sex
lr  <- as.numeric(2 * (logLik(fit_joint) - logLik(fit_base)))
dfd <- attr(logLik(fit_joint), "df") - attr(logLik(fit_base), "df")
lr_p <- pchisq(lr, df = dfd, lower.tail = FALSE)

fit_null <- glm(Y01 ~ 1, family = binomial, data = d)
ll_joint <- as.numeric(logLik(fit_joint))
ll_null  <- as.numeric(logLik(fit_null))
n <- nrow(d)
mcfadden <- 1 - ll_joint / ll_null
cox      <- 1 - exp((2 / n) * (ll_null - ll_joint))
nag      <- cox / (1 - exp((2 / n) * ll_null))

co <- summary(fit_joint)$coefficients
per_prs_rows <- lapply(BASES, function(base) {
    term <- sprintf("PRS_%s_z", base)
    if (!(term %in% rownames(co))) return(NULL)
    est <- co[term, "Estimate"]; se <- co[term, "Std. Error"]
    z   <- co[term, "z value"];  P  <- co[term, "Pr(>|z|)"]
    data.frame(PRS = base,
               beta_z = signif(est, 4),
               SE     = signif(se, 4),
               OR_perSD  = signif(exp(est), 4),
               OR_CI_lo  = signif(exp(est - 1.96 * se), 4),
               OR_CI_hi  = signif(exp(est + 1.96 * se), 4),
               z         = signif(z, 4),
               P         = signif(P, 4),
               best_Pt   = best_pt_by_base[base],
               N_SNP     = nsnp_by_base[base],
               PRSice_R2 = signif(r2_by_base[base], 4),
               stringsAsFactors = FALSE)
})
coef_tab <- do.call(rbind, per_prs_rows)

model_tab <- data.frame(
    N              = n,
    N_cases        = sum(d$Y01 == 1),
    N_ctrls        = sum(d$Y01 == 0),
    LR_chisq_PRS   = signif(lr,  4),
    LR_df          = dfd,
    LR_P           = signif(lr_p, 4),
    McFadden_R2    = signif(mcfadden, 4),
    CoxSnell_R2    = signif(cox,      4),
    Nagelkerke_R2  = signif(nag,      4),
    AIC_joint      = signif(AIC(fit_joint), 6),
    AIC_baseline   = signif(AIC(fit_base),  6),
    stringsAsFactors = FALSE)

# AUC via Mann-Whitney U, plus the ROC trajectory
auc_and_roc <- function(y, pred) {
    p1 <- pred[y == 1]; p0 <- pred[y == 0]
    n1 <- length(p1);   n0 <- length(p0)
    r  <- rank(c(p1, p0))
    u  <- sum(r[seq_len(n1)]) - n1 * (n1 + 1) / 2
    auc <- u / (n1 * n0)
    o   <- order(pred, decreasing = TRUE)
    ys  <- y[o]
    tpr <- c(0, cumsum(ys == 1) / n1)
    fpr <- c(0, cumsum(ys == 0) / n0)
    list(auc = auc, tpr = tpr, fpr = fpr)
}
pred_joint    <- predict(fit_joint,    type = "response")
pred_base     <- predict(fit_base,     type = "response")
pred_base_pcs <- predict(fit_base_pcs, type = "response")
pred_base_sex <- predict(fit_base_sex, type = "response")
roc_joint    <- auc_and_roc(d$Y01, pred_joint)
roc_base     <- auc_and_roc(d$Y01, pred_base)
roc_base_pcs <- auc_and_roc(d$Y01, pred_base_pcs)
roc_base_sex <- auc_and_roc(d$Y01, pred_base_sex)
model_tab$AUC_joint       <- signif(roc_joint$auc,    4)
model_tab$AUC_base_PCsSex <- signif(roc_base$auc,     4)
model_tab$AUC_base_PCs    <- signif(roc_base_pcs$auc, 4)
model_tab$AUC_base_Sex    <- signif(roc_base_sex$auc, 4)
model_tab$AUC_gain        <- signif(roc_joint$auc - roc_base$auc, 4)

coef_out  <- file.path(WORK, paste0(PFX, "_coef.tsv"))
model_out <- file.path(WORK, paste0(PFX, "_model.tsv"))
write.table(coef_tab,  coef_out,  sep = "\t", row.names = FALSE, quote = FALSE)
write.table(model_tab, model_out, sep = "\t", row.names = FALSE, quote = FALSE)
cat("\nPer-PRS coefficients (joint model):\n"); print(coef_tab,  row.names = FALSE)
cat("\nModel-level stats:\n");                  print(model_tab, row.names = FALSE)
cat("\nWrote ", coef_out,  "\n", sep = "")
cat("Wrote ",   model_out, "\n", sep = "")

lab_joint  <- sprintf("Joint %d-PRS + PCs + Sex (AUC=%.3f)", NB, roc_joint$auc)
lab_base   <- sprintf("PCs + Sex (AUC=%.3f)",               roc_base$auc)
lab_pcs    <- sprintf("PCs only  (AUC=%.3f)",               roc_base_pcs$auc)
lab_sex    <- sprintf("Sex only  (AUC=%.3f)",               roc_base_sex$auc)
roc_df <- rbind(
    data.frame(fpr = roc_joint$fpr,    tpr = roc_joint$tpr,    Model = lab_joint),
    data.frame(fpr = roc_base$fpr,     tpr = roc_base$tpr,     Model = lab_base),
    data.frame(fpr = roc_base_pcs$fpr, tpr = roc_base_pcs$tpr, Model = lab_pcs),
    data.frame(fpr = roc_base_sex$fpr, tpr = roc_base_sex$tpr, Model = lab_sex)
)
roc_df$Model <- factor(roc_df$Model,
                       levels = c(lab_joint, lab_base, lab_pcs, lab_sex))

g_roc <- ggplot(roc_df, aes(fpr, tpr, color = Model, linetype = Model)) +
    geom_abline(intercept = 0, slope = 1, color = "grey70", linetype = "dotted") +
    geom_line(linewidth = 0.9) +
    scale_color_manual(values = setNames(c("#332288","grey30","#CC6677","#3E8E7E"),
                                          levels(roc_df$Model))) +
    scale_linetype_manual(values = setNames(c("solid","dashed","dashed","dashed"),
                                            levels(roc_df$Model))) +
    coord_equal() +
    labs(title = sprintf("Within-case %s: joint %d-PRS model\nwith diagnostic baselines", OUTCOME, NB),
         subtitle = sprintf(
             "LR test (%d PRSes over PCs+Sex): chisq=%.2f (df=%d), P=%s\nNagelkerke R2=%.3f",
             NB, lr, dfd, signif(lr_p, 3), nag),
         x = "False positive rate",
         y = "True positive rate") +
    theme_minimal(base_size = 12) +
    theme(plot.title    = element_text(face = "bold", size = 13,
                                       margin = margin(b = 4)),
          plot.subtitle = element_text(size = 10, color = "grey30",
                                       margin = margin(b = 6)),
          plot.title.position = "plot",
          plot.margin   = margin(t = 8, r = 12, b = 6, l = 8),
          legend.position = "bottom",
          legend.title    = element_blank(),
          legend.text     = element_text(size = 9),
          panel.grid.minor = element_blank()) +
    guides(color = guide_legend(ncol = 2), linetype = guide_legend(ncol = 2))
ggsave(file.path(WORK, paste0(PFX, "_ROC.png")), g_roc,
       width = 8.5, height = 7.5, dpi = 300)
ggsave(file.path(WORK, paste0(PFX, "_ROC.pdf")), g_roc,
       width = 8.5, height = 7.5)
cat("Wrote ", PFX, "_ROC.{png,pdf}\n", sep = "")

fp <- coef_tab
fp$PRS <- factor(tools::toTitleCase(fp$PRS),
                 levels = tools::toTitleCase(BASES))
fp$Sig <- ifelse(fp$P < 0.001, "***",
          ifelse(fp$P < 0.01,  "**",
          ifelse(fp$P < 0.05,  "*",
          ifelse(fp$P < 0.1,   ".", ""))))

g_forest <- ggplot(fp, aes(x = OR_perSD, y = PRS)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
    geom_errorbarh(aes(xmin = OR_CI_lo, xmax = OR_CI_hi),
                   height = 0.22, linewidth = 0.6, color = "grey30") +
    geom_point(size = 3, color = "#332288") +
    geom_text(aes(label = sprintf("OR=%.2f [%.2f, %.2f]   P=%.3g %s",
                                   OR_perSD, OR_CI_lo, OR_CI_hi, P, Sig),
                  x = OR_CI_hi),
              hjust = -0.1, size = 3.2, color = "grey15") +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.55))) +
    labs(title = sprintf("Within-case %s: per-SD OR from the joint %d-PRS model", OUTCOME, NB),
         subtitle = "Adjusted for PC1-4 and Sex. All PRSes z-scored within the cases.",
         x = "Odds ratio per SD of PRS (joint model)",
         y = NULL) +
    theme_minimal(base_size = 12) +
    theme(plot.title    = element_text(face = "bold"),
          plot.subtitle = element_text(size = 10, color = "grey30"),
          panel.grid.minor = element_blank())
ggsave(file.path(WORK, paste0(PFX, "_forest.png")), g_forest,
       width = 8.5, height = 4.5, dpi = 300)
ggsave(file.path(WORK, paste0(PFX, "_forest.pdf")), g_forest,
       width = 8.5, height = 4.5)
cat("Wrote ", PFX, "_forest.{png,pdf}\n", sep = "")
