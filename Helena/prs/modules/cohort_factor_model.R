#!/usr/bin/env Rscript
# Contrasts of cases vs each control cohort from one model per score: PRS_z ~ Cohort + PC1-4 + Sex.
# Usage: Rscript cohort_factor_model.R [WORKDIR]   (reads results/ from 06_run_prsice.sh)
suppressPackageStartupMessages({ library(ggplot2) })

args       <- commandArgs(trailingOnly = TRUE)
WORK_DIR   <- if (length(args) >= 1) args[1] else Sys.getenv("WORKDIR", "/re_gecip/<your_gecip>/<project>/<cohort>")
PRS_LIST   <- strsplit(Sys.getenv("PRS_LIST", "epilepsy intelligence mvage"), "[[:space:]]+")[[1]]
CASE_LABEL <- Sys.getenv("CASE_LABEL", "Cases")
RES_DIR    <- file.path(WORK_DIR, "results")

PHENOS <- PRS_LIST[file.exists(file.path(RES_DIR, sprintf("PRS_%s.summary", PRS_LIST)))]
cat("Phenotypes:", paste(PHENOS, collapse = ", "), "\n")

setwd(WORK_DIR)

read_ids <- function(p) {
  if (!file.exists(p)) return(data.frame(FID = character(), IID = character()))
  read.table(p, header = FALSE, col.names = c("FID","IID"),
             stringsAsFactors = FALSE)
}
cs  <- read_ids("merged/case_samples.txt")
dee <- read_ids("merged/dee_samples.txt")
epi <- read_ids("merged/epi_samples.txt")
hc  <- read_ids("merged/hc_samples.txt")
cov <- read.table("cov.txt", header = TRUE, stringsAsFactors = FALSE)

label_cohort <- function(iid) {
  ifelse(iid %in% cs$IID,  CASE_LABEL,
  ifelse(iid %in% dee$IID, "DEE",
  ifelse(iid %in% epi$IID, "Epi",
  ifelse(iid %in% hc$IID,  "HC", NA_character_))))
}

load_prs <- function(p) {
  scores_file <- file.path(RES_DIR, sprintf("PRS_%s.all.score", p))
  summ_file   <- file.path(RES_DIR, sprintf("PRS_%s.summary",   p))
  if (!file.exists(scores_file)) { message("missing ", scores_file); return(NULL) }

  scores <- read.table(scores_file, header = TRUE, check.names = FALSE,
                       stringsAsFactors = FALSE)
  if (!"IID" %in% names(scores)) names(scores)[1:2] <- c("FID","IID")
  pt_cols <- setdiff(names(scores), c("FID","IID"))
  pt_nums <- suppressWarnings(as.numeric(sub("^Pt_", "", pt_cols)))

  best_pt <- NA
  if (file.exists(summ_file)) {
    s <- tryCatch(read.table(summ_file, header = TRUE), error = function(e) NULL)
    if (!is.null(s)) best_pt <- s$Threshold[1]
  }
  best_col <- NA_character_
  if (!is.na(best_pt) && any(!is.na(pt_nums))) {
    idx <- which.min(abs(pt_nums - best_pt))
    if (length(idx) == 1 && !is.na(pt_nums[idx]) &&
        abs(pt_nums[idx] - best_pt) < 1e-9 + 1e-3 * abs(best_pt))
      best_col <- pt_cols[idx]
  }
  if (is.na(best_col)) {
    best_file <- file.path(RES_DIR, sprintf("PRS_%s.best", p))
    if (file.exists(best_file)) {
      bf <- read.table(best_file, header = TRUE, stringsAsFactors = FALSE)
      scores <- data.frame(FID = bf$FID, IID = bf$IID, PRS_best = bf$PRS,
                           stringsAsFactors = FALSE)
      best_col <- "PRS_best"
    } else if (length(pt_cols)) {
      best_col <- pt_cols[length(pt_cols)]
    } else {
      return(NULL)
    }
  }
  cat(sprintf("[%s] best Pt=%s  col=%s\n", p, best_pt, best_col))

  d <- data.frame(IID = scores$IID, PRS = scores[[best_col]],
                  stringsAsFactors = FALSE)
  d <- merge(d, cov, by = "IID", all.x = TRUE)
  d$Cohort <- label_cohort(d$IID)
  d <- d[!is.na(d$Cohort) & !is.na(d$PRS), ]
  d$PRS_z <- as.numeric(scale(d$PRS))
  d
}

CONTRASTS <- paste(CASE_LABEL, "vs", c("GEL Epilepsy controls", "GEL Healthy controls", "DEE cases"))

# Epi is the reference; refitting with HC and DEE as reference gives the other two contrasts
get_contrasts <- function(p) {
  d <- load_prs(p)
  if (is.null(d) || nrow(d) == 0) return(NULL)

  d$Cohort <- factor(d$Cohort, levels = c("Epi","HC","DEE",CASE_LABEL))
  case_term <- paste0("Cohort", CASE_LABEL)

  pull <- function(fit, term) {
    s <- summary(fit)$coefficients
    if (!(term %in% rownames(s))) return(c(NA, NA, NA, NA))
    s[term, c("Estimate","Std. Error","t value","Pr(>|t|)")]
  }

  fit_epi <- lm(PRS_z ~ Cohort + PC1 + PC2 + PC3 + PC4 + factor(Sex), data = d)
  vs_epi  <- pull(fit_epi, case_term)

  d2 <- d; d2$Cohort <- relevel(d2$Cohort, ref = "HC")
  fit_hc <- lm(PRS_z ~ Cohort + PC1 + PC2 + PC3 + PC4 + factor(Sex), data = d2)
  vs_hc  <- pull(fit_hc, case_term)

  d3 <- d; d3$Cohort <- relevel(d3$Cohort, ref = "DEE")
  fit_dee <- lm(PRS_z ~ Cohort + PC1 + PC2 + PC3 + PC4 + factor(Sex), data = d3)
  vs_dee  <- pull(fit_dee, case_term)

  data.frame(
    Phenotype = p,
    Contrast  = CONTRASTS,
    Effect    = c(vs_epi[1], vs_hc[1], vs_dee[1]),
    SE        = c(vs_epi[2], vs_hc[2], vs_dee[2]),
    t         = c(vs_epi[3], vs_hc[3], vs_dee[3]),
    P         = c(vs_epi[4], vs_hc[4], vs_dee[4]),
    N         = nrow(d),
    stringsAsFactors = FALSE
  )
}

out <- do.call(rbind, lapply(PHENOS, get_contrasts))
out$CI_lo <- out$Effect - 1.96 * out$SE
out$CI_hi <- out$Effect + 1.96 * out$SE
out <- out[, c("Phenotype","Contrast","Effect","SE","CI_lo","CI_hi","t","P","N")]

write.table(out, file.path(RES_DIR, "cohort_contrasts_table.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("\nWrote cohort_contrasts_table.tsv\n")
print(out)

wrap_label <- function(x, width = 20) {
  vapply(as.character(x), function(s) {
    if (nchar(s) <= width) return(s)
    spaces <- gregexpr(" ", s)[[1]]
    if (length(spaces) == 0 || spaces[1] < 0) return(s)
    mid  <- nchar(s) / 2
    best <- spaces[which.min(abs(spaces - mid))]
    paste0(substr(s, 1, best - 1), "\n", substr(s, best + 1, nchar(s)))
  }, character(1))
}

out$Phenotype <- factor(out$Phenotype, levels = unique(out$Phenotype))
out$Contrast  <- factor(out$Contrast, levels = CONTRASTS)
out$Sig <- ifelse(out$P < 0.001, "***",
           ifelse(out$P < 0.01,  "**",
           ifelse(out$P < 0.05,  "*",
           ifelse(out$P < 0.1,   ".", ""))))

g <- ggplot(out, aes(x = Effect, y = Contrast, color = Contrast)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbarh(aes(xmin = CI_lo, xmax = CI_hi), height = 0.25, linewidth = 0.6) +
  geom_point(size = 2.5) +
  geom_text(aes(label = Sig, x = CI_hi), hjust = -0.3, size = 4, color = "grey15") +
  facet_wrap(~ Phenotype, scales = "free_x", ncol = 1) +
  scale_color_manual(values = setNames(c("#E69F00", "#3498DB", "#9B59B6"), CONTRASTS),
                     guide = "none") +
  scale_y_discrete(labels = wrap_label) +
  labs(title = sprintf("%s vs each control cohort, by phenotype", CASE_LABEL),
       subtitle = "Cohort-factor model: PRS_z ~ Cohort + PC1-4 + Sex (4 cohorts)",
       x = sprintf("Effect on standardised PRS (%s - control)", CASE_LABEL),
       y = NULL) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"))

ggsave(file.path(RES_DIR, "cohort_contrasts_forest.png"), g,
       width = 8, height = 1.6 + 1.6 * length(unique(out$Phenotype)), dpi = 300)
ggsave(file.path(RES_DIR, "cohort_contrasts_forest.pdf"), g,
       width = 8, height = 1.6 + 1.6 * length(unique(out$Phenotype)))
cat("Wrote cohort_contrasts_forest.{png,pdf}\n")
