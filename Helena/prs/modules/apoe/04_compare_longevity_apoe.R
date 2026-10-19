#!/usr/bin/env Rscript
# Pairwise cohort P values for the longevity PRS with and without the APOE region, side by side.
# Usage: Rscript 04_compare_longevity_apoe.R [WORKDIR]   (needs 03_prsice_longevity_no_apoe.sh)
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })

args       <- commandArgs(trailingOnly = TRUE)
WORK_DIR   <- if (length(args) >= 1) args[1] else Sys.getenv("WORKDIR", "/re_gecip/<your_gecip>/<project>/<cohort>")
SCORE      <- Sys.getenv("LONGEVITY_SCORE", "longevity")
CASE_LABEL <- Sys.getenv("CASE_LABEL", "Cases")
RES_DIR    <- file.path(WORK_DIR, "results")
OUT_DIR    <- file.path(WORK_DIR, "apoe")

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
cov <- read.table(file.path(WORK_DIR, "cov.txt"),
                  header=TRUE, stringsAsFactors=FALSE)

load_prs <- function(score_f, summ_f) {
  if (!file.exists(score_f)) stop("missing ", score_f)
  x <- read.table(score_f, header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)
  if (!"IID" %in% names(x)) names(x)[1:2] <- c("FID","IID")
  pt_cols <- setdiff(names(x), c("FID","IID"))
  pt_nums <- suppressWarnings(as.numeric(sub("^Pt_", "", pt_cols)))
  best_pt <- NA_real_
  if (file.exists(summ_f)) {
    s <- tryCatch(read.table(summ_f, header=TRUE), error=function(e) NULL)
    if (!is.null(s)) best_pt <- s$Threshold[1]
  }
  best_col <- if (!is.na(best_pt) && any(!is.na(pt_nums))) {
    i <- which.min(abs(pt_nums - best_pt)); pt_cols[i]
  } else pt_cols[length(pt_cols)]
  cat(sprintf("  %s : Pt=%s -> %s\n", basename(score_f), best_pt, best_col))
  data.frame(IID=x$IID, PRS=x[[best_col]], stringsAsFactors=FALSE)
}

prs_with <- load_prs(file.path(RES_DIR, sprintf("PRS_%s.all.score", SCORE)),
                     file.path(RES_DIR, sprintf("PRS_%s.summary", SCORE)))
prs_no   <- load_prs(file.path(OUT_DIR, sprintf("PRS_%s_noAPOE.all.score", SCORE)),
                     file.path(OUT_DIR, sprintf("PRS_%s_noAPOE.summary", SCORE)))

prep <- function(prs) {
  x <- merge(prs, cov, by="IID", all.x=TRUE)
  x$Cohort <- label_cohort(x$IID)
  x <- x[!is.na(x$Cohort) & !is.na(x$PRS), ]
  x$PRS_z <- as.numeric(scale(x$PRS))
  fit <- lm(PRS ~ PC1+PC2+PC3+PC4+factor(Sex), data=x)
  x$PRS_adj <- as.numeric(scale(residuals(fit)))
  x
}
d_with <- prep(prs_with)
d_no   <- prep(prs_no)

lvls <- c(CASE_LABEL,"GEL Epilepsy controls","DEE cases","GEL Healthy controls")
lvls <- lvls[lvls %in% intersect(unique(d_with$Cohort), unique(d_no$Cohort))]

pairwise_p <- function(dat, ref, other) {
  dat$Cohort <- relevel(factor(dat$Cohort, levels=lvls), ref=ref)
  fit <- tryCatch(lm(PRS_z ~ Cohort + PC1+PC2+PC3+PC4+factor(Sex), data=dat),
                  error=function(e) NULL)
  if (is.null(fit)) return(NA_real_)
  s <- summary(fit)$coefficients
  term <- paste0("Cohort", other)
  if (!term %in% rownames(s)) return(NA_real_)
  s[term, "Pr(>|t|)"]
}
pairs <- combn(lvls, 2, simplify=FALSE)
res <- data.frame(
  cohort_A = sapply(pairs, `[`, 1),
  cohort_B = sapply(pairs, `[`, 2),
  p_with_APOE = sapply(pairs, function(p) pairwise_p(d_with, p[1], p[2])),
  p_no_APOE   = sapply(pairs, function(p) pairwise_p(d_no,   p[1], p[2]))
)
res$conclusion_unchanged <- with(res, (p_with_APOE < 0.05) == (p_no_APOE < 0.05))
write.table(res, file.path(OUT_DIR, sprintf("PRS_%s_APOE_sensitivity.tsv", SCORE)),
            sep="\t", row.names=FALSE, quote=FALSE)
cat("\nPairwise APOE sensitivity:\n"); print(res, row.names=FALSE)

cohort_cols <- setNames(c("#3498DB", "#E69F00", "#9B59B6", "#332288"),
                        c("GEL Healthy controls", "GEL Epilepsy controls", "DEE cases", CASE_LABEL))

fmt_p <- function(x) {
  if (is.na(x)) return("NA")
  if (x < 1e-4) return(format(x, scientific=TRUE, digits=2))
  formatC(x, format="fg", digits=2)
}

# same y-limit in both panels so the brackets sit at comparable heights
y_top <- max(c(d_with$PRS_adj, d_no$PRS_adj), na.rm=TRUE)
step  <- 0.45
n_pairs <- length(pairs)
y_max <- y_top + step * (n_pairs + 1)

mk_plot <- function(dat, ttl) {
  dat$Cohort <- factor(dat$Cohort, levels=lvls)
  ns <- table(dat$Cohort)
  labs <- setNames(sprintf("%s\n(n=%d)", names(ns), as.integer(ns)), names(ns))
  ps <- sapply(pairs, function(p) pairwise_p(dat, p[1], p[2]))
  br <- data.frame(
    x1 = sapply(pairs, function(p) match(p[1], lvls)),
    x2 = sapply(pairs, function(p) match(p[2], lvls)),
    p  = ps
  )
  br <- br[order(br$x2 - br$x1), ]
  br$y <- y_top + step * seq_len(nrow(br))
  br$label <- vapply(br$p, fmt_p, character(1))
  ggplot(dat, aes(Cohort, PRS_adj, color=Cohort)) +
    geom_hline(yintercept=0, linetype="dashed", color="grey40") +
    geom_jitter(width=0.28, alpha=0.5, size=0.8) +
    geom_boxplot(fill=NA, outlier.shape=NA, width=0.55, linewidth=0.6) +
    geom_segment(data=br, inherit.aes=FALSE,
                 aes(x=x1, xend=x2, y=y, yend=y), color="grey25") +
    geom_segment(data=br, inherit.aes=FALSE,
                 aes(x=x1, xend=x1, y=y, yend=y-0.1), color="grey25") +
    geom_segment(data=br, inherit.aes=FALSE,
                 aes(x=x2, xend=x2, y=y, yend=y-0.1), color="grey25") +
    geom_text(data=br, inherit.aes=FALSE,
              aes(x=(x1+x2)/2, y=y+0.12, label=label),
              size=3.0, color="grey15") +
    scale_color_manual(values=cohort_cols, guide="none") +
    scale_x_discrete(labels=labs) +
    coord_cartesian(ylim=c(NA, y_max)) +
    labs(title=ttl, x=NULL, y=sprintf("PRS - %s (adj, z)", tools::toTitleCase(SCORE))) +
    theme_minimal(base_size=11)
}
g <- mk_plot(d_with, "With APOE") + mk_plot(d_no, "Without APOE (sensitivity)")
ggsave(file.path(OUT_DIR, sprintf("PRS_%s_APOE_sensitivity.png", SCORE)),
       g, width=12, height=6.5, dpi=300)
cat("Done. Outputs in ", OUT_DIR, "\n", sep="")
