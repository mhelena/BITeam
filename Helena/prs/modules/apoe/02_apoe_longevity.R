#!/usr/bin/env Rscript
# APOE genotype by cohort, and the longevity PRS stratified by and adjusted for APOE.
# Usage: Rscript 02_apoe_longevity.R [WORKDIR]   (needs 01_extract_apoe_status.sh and the longevity PRSice run)
suppressPackageStartupMessages({ library(ggplot2) })

args       <- commandArgs(trailingOnly = TRUE)
WORK_DIR   <- if (length(args) >= 1) args[1] else Sys.getenv("WORKDIR", "/re_gecip/<your_gecip>/<project>/<cohort>")
SCORE      <- Sys.getenv("LONGEVITY_SCORE", "longevity")
CASE_LABEL <- Sys.getenv("CASE_LABEL", "Cases")
RES_DIR    <- file.path(WORK_DIR, "results")
OUT_DIR    <- file.path(WORK_DIR, "apoe")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

read_ids <- function(p) {
  if (!file.exists(p)) return(data.frame(FID=character(), IID=character()))
  read.table(p, header=FALSE, col.names=c("FID","IID"), stringsAsFactors=FALSE)
}
cs  <- read_ids(file.path(WORK_DIR, "merged/case_samples.txt"))
dee <- read_ids(file.path(WORK_DIR, "merged/dee_samples.txt"))
epi <- read_ids(file.path(WORK_DIR, "merged/epi_samples.txt"))
hc  <- read_ids(file.path(WORK_DIR, "merged/hc_samples.txt"))
label_cohort <- function(iid) {
  ifelse(iid %in% cs$IID,  CASE_LABEL,
  ifelse(iid %in% dee$IID, "DEE cases",
  ifelse(iid %in% epi$IID, "GEL Epilepsy controls",
  ifelse(iid %in% hc$IID,  "GEL Healthy controls", NA_character_))))
}

apoe <- read.table(file.path(OUT_DIR, "APOE_status.tsv"),
                   header=TRUE, sep="\t", stringsAsFactors=FALSE, na.strings="NA")
cov <- read.table(file.path(WORK_DIR, "cov.txt"),
                  header=TRUE, stringsAsFactors=FALSE)

score_f <- file.path(RES_DIR, sprintf("PRS_%s.all.score", SCORE))
summ_f  <- file.path(RES_DIR, sprintf("PRS_%s.summary", SCORE))
if (!file.exists(score_f)) stop("missing ", score_f)
scores <- read.table(score_f, header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)
if (!"IID" %in% names(scores)) names(scores)[1:2] <- c("FID","IID")
pt_cols <- setdiff(names(scores), c("FID","IID"))
pt_nums <- suppressWarnings(as.numeric(sub("^Pt_", "", pt_cols)))
best_pt <- NA_real_
if (file.exists(summ_f)) {
  s <- tryCatch(read.table(summ_f, header=TRUE), error=function(e) NULL)
  if (!is.null(s)) best_pt <- s$Threshold[1]
}
best_col <- if (!is.na(best_pt) && any(!is.na(pt_nums))) {
  i <- which.min(abs(pt_nums - best_pt)); pt_cols[i]
} else pt_cols[length(pt_cols)]
cat(sprintf("Using PRS threshold Pt = %s (column %s)\n", best_pt, best_col))

d <- data.frame(IID=scores$IID, PRS=scores[[best_col]], stringsAsFactors=FALSE)
d <- merge(d, cov, by="IID", all.x=TRUE)
d <- merge(d, apoe[, c("IID","APOE_genotype","APOE_e4_dose","APOE_e2_dose")],
           by="IID", all.x=TRUE)
d$Cohort <- label_cohort(d$IID)
d <- d[!is.na(d$Cohort) & !is.na(d$PRS), ]

d_all <- d
d <- d[!is.na(d$APOE_genotype) & d$APOE_genotype != "excluded_double_het", ]
n_excl <- nrow(d_all) - nrow(d)
cat(sprintf("Excluded %d samples with ambiguous double-het (APOE)\n", n_excl))

lvls <- c(CASE_LABEL,"GEL Epilepsy controls","DEE cases","GEL Healthy controls")
lvls <- lvls[lvls %in% unique(d$Cohort)]
d$Cohort <- factor(d$Cohort, levels=lvls)

# APOE genotype x cohort
tab <- table(d$APOE_genotype, d$Cohort)
tab <- tab[intersect(c("e2/e2","e2/e3","e3/e3","e3/e4","e4/e4"), rownames(tab)),
           , drop=FALSE]
pct <- sweep(tab, 2, colSums(tab), FUN="/") * 100
tab_out <- as.data.frame.matrix(tab)
for (co in colnames(pct)) tab_out[[co]] <- sprintf("%d (%.1f%%)", tab[, co], pct[, co])
tab_out <- cbind(APOE_genotype=rownames(tab_out), tab_out)
chi <- suppressWarnings(chisq.test(tab))
cat(sprintf("Chi-square across genotypes x cohorts: X2=%.2f, df=%d, p=%.3g\n",
            chi$statistic, chi$parameter, chi$p.value))
write.table(tab_out, file.path(OUT_DIR, "Table1_APOE_by_cohort.tsv"),
            sep="\t", row.names=FALSE, quote=FALSE)

# Within each cohort: PRS in e4 carriers vs e2 carriers
ttest_rows <- list()
for (co in lvls) {
  sub <- d[d$Cohort == co, ]
  e4 <- sub$PRS[sub$APOE_genotype %in% c("e3/e4","e4/e4")]
  e2 <- sub$PRS[sub$APOE_genotype %in% c("e2/e2","e2/e3")]
  if (length(e4) >= 2 && length(e2) >= 2) {
    tt <- t.test(e4, e2)
    ttest_rows[[co]] <- data.frame(
      Cohort=co, n_e4=length(e4), n_e2=length(e2),
      mean_PRS_e4=mean(e4), mean_PRS_e2=mean(e2),
      t=unname(tt$statistic), df=unname(tt$parameter), p=tt$p.value
    )
  } else {
    ttest_rows[[co]] <- data.frame(Cohort=co, n_e4=length(e4), n_e2=length(e2),
      mean_PRS_e4=mean(e4), mean_PRS_e2=mean(e2), t=NA, df=NA, p=NA)
  }
}
ttest_df <- do.call(rbind, ttest_rows)
write.table(ttest_df, file.path(OUT_DIR, "APOE_carrier_ttests.tsv"),
            sep="\t", row.names=FALSE, quote=FALSE)

# Pairwise cohort P values, with and without e4/e2 dose in the model
d$e4_dose <- as.integer(d$APOE_e4_dose)
d$e2_dose <- as.integer(d$APOE_e2_dose)
d$PRS_z <- as.numeric(scale(d$PRS))

pairwise_p <- function(ref, other, adj_apoe) {
  d2 <- d; d2$Cohort <- relevel(d2$Cohort, ref=ref)
  f <- if (adj_apoe)
    PRS_z ~ Cohort + e4_dose + e2_dose + PC1+PC2+PC3+PC4+factor(Sex)
  else
    PRS_z ~ Cohort + PC1+PC2+PC3+PC4+factor(Sex)
  fit <- tryCatch(lm(f, data=d2), error=function(e) NULL)
  if (is.null(fit)) return(NA_real_)
  s <- summary(fit)$coefficients
  term <- paste0("Cohort", other)
  if (!term %in% rownames(s)) return(NA_real_)
  s[term, "Pr(>|t|)"]
}
pairs <- combn(lvls, 2, simplify=FALSE)
lm_rows <- data.frame(
  cohort_A = sapply(pairs, `[`, 1),
  cohort_B = sapply(pairs, `[`, 2),
  p_unadjusted = sapply(pairs, function(p) pairwise_p(p[1], p[2], adj_apoe=FALSE)),
  p_APOE_adjusted = sapply(pairs, function(p) pairwise_p(p[1], p[2], adj_apoe=TRUE))
)
write.table(lm_rows, file.path(OUT_DIR, "APOE_adjusted_lm.tsv"),
            sep="\t", row.names=FALSE, quote=FALSE)

d$APOE_genotype <- factor(d$APOE_genotype,
                          levels=intersect(c("e2/e2","e2/e3","e3/e3","e3/e4","e4/e4"),
                                           unique(d$APOE_genotype)))
geno_lvls <- levels(d$APOE_genotype)
cohort_cols <- setNames(c("#3498DB", "#E69F00", "#9B59B6", "#332288"),
                        c("GEL Healthy controls", "GEL Epilepsy controls", "DEE cases", CASE_LABEL))
fit <- lm(PRS ~ PC1+PC2+PC3+PC4+factor(Sex), data=d)
d$PRS_adj <- as.numeric(scale(residuals(fit)))

fmt_p <- function(x) {
  if (is.na(x))     return("NA")
  if (x < 1e-4)     return(format(x, scientific=TRUE, digits=2))
  formatC(x, format="fg", digits=2)
}
# cases vs one control within an APOE stratum, adjusted for PCs and sex
case_vs_p <- function(sub, other) {
  n_case <- sum(sub$Cohort==CASE_LABEL); n_oth <- sum(sub$Cohort==other)
  if (n_case < 3 || n_oth < 3) return(NA_real_)
  sub$Cohort <- factor(sub$Cohort, levels=lvls)
  sub$Cohort <- relevel(sub$Cohort, ref=CASE_LABEL)
  fit <- tryCatch(lm(PRS_z ~ Cohort + PC1+PC2+PC3+PC4+factor(Sex), data=sub),
                  error=function(e) NULL)
  if (is.null(fit)) return(NA_real_)
  s <- summary(fit)$coefficients
  term <- paste0("Cohort", other)
  if (!term %in% rownames(s)) return(NA_real_)
  s[term, "Pr(>|t|)"]
}
p_lab <- do.call(rbind, lapply(geno_lvls, function(g0) {
  sub <- d[d$APOE_genotype == g0, ]
  p_epi <- case_vs_p(sub, "GEL Epilepsy controls")
  p_dee <- case_vs_p(sub, "DEE cases")
  p_hc  <- case_vs_p(sub, "GEL Healthy controls")
  label <- paste(
    sprintf("vs Epi: p=%s", fmt_p(p_epi)),
    sprintf("vs DEE: p=%s", fmt_p(p_dee)),
    sprintf("vs HC:  p=%s", fmt_p(p_hc)),
    sep="\n"
  )
  data.frame(APOE_genotype=g0, label=label, stringsAsFactors=FALSE)
}))
p_lab$APOE_genotype <- factor(p_lab$APOE_genotype, levels=geno_lvls)
write.table(p_lab, file.path(OUT_DIR, "APOE_boxplot_pvalues.tsv"),
            sep="\t", row.names=FALSE, quote=FALSE)

y_top <- max(d$PRS_adj, na.rm=TRUE)
y_lab <- y_top + 0.9
g <- ggplot(d, aes(APOE_genotype, PRS_adj, color=Cohort)) +
  geom_hline(yintercept=0, linetype="dashed", color="grey40") +
  geom_boxplot(outlier.shape=NA, position=position_dodge(width=0.75), width=0.6) +
  geom_jitter(position=position_jitterdodge(jitter.width=0.15, dodge.width=0.75),
              alpha=0.5, size=0.8) +
  geom_text(data=p_lab, inherit.aes=FALSE,
            aes(x=APOE_genotype, y=y_lab, label=label),
            size=2.6, color="grey20", lineheight=0.95, vjust=0) +
  scale_color_manual(values=cohort_cols) +
  coord_cartesian(ylim=c(NA, y_lab + 1.0), clip="off") +
  labs(title=sprintf("%s PRS by APOE genotype and cohort", tools::toTitleCase(SCORE)),
       subtitle=sprintf("P-values above each stratum: %s vs each control (lm adj. PC1-4 + Sex)", CASE_LABEL),
       x="APOE genotype",
       y=sprintf("PRS - %s (adjusted, standardised)", tools::toTitleCase(SCORE))) +
  theme_minimal(base_size=12) +
  theme(plot.subtitle=element_text(size=9, color="grey30"))
ggsave(file.path(OUT_DIR, "APOE_boxplot_by_cohort.png"),
       g, width=10, height=6, dpi=300)

cat("Done. Outputs in ", OUT_DIR, "\n", sep="")
