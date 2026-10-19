#!/usr/bin/env Rscript
# PRS tertile/quartile stratification of pBAG and PRS x clinical factor interactions.
# Usage: Rscript scripts/06_M5_stratification.R config.sh (needs PRS_pBAG_linear_results.txt from step 4)

suppressPackageStartupMessages({
  library(ggplot2); library(ggpubr); library(dplyr); library(tidyr); library(reshape2)
})

args <- commandArgs(trailingOnly = FALSE)
source(file.path(dirname(sub("^--file=", "", grep("^--file=", args, value = TRUE))), "config.R"))
cfg <- read_config()

META <- cfg("METADATA")
COV  <- cfg("COVARIATES")
SCO  <- cfg("PRS_DIR")
OUT  <- cfg("OUT_DIR")
M3   <- file.path(OUT, "PRS_pBAG_linear_results.txt")

meta <- read.delim(META, stringsAsFactors = FALSE, check.names = FALSE)
meta <- meta[!(meta$MRI_findings %in% c("Post surgery", "Large lesion")), ]
cat("After exclusion: n =", nrow(meta), "\n")

names(meta)[names(meta) == "Age at onset"]            <- "AgeOnset"
names(meta)[names(meta) == "EpilepsySyndrome_forPRS"] <- "Syndrome"
names(meta)[names(meta) == "HistoryOfConvulsions"]    <- "TC_seizures"
meta$Syndrome    <- ifelse(meta$Syndrome %in% c("Focal","Generalised"), meta$Syndrome, NA)
meta$TC_seizures <- ifelse(meta$TC_seizures == "History of convulsions", 1,
                    ifelse(meta$TC_seizures == "NO History of convulsions", 0, NA))
meta$DRE_bin     <- ifelse(meta$DRE == "DRE at time of scan", 1,
                    ifelse(meta$DRE == "NO DRE at time of scan", 0, NA))
meta$VPA_bin     <- ifelse(meta$VPA == "VPA exposure", 1,
                    ifelse(meta$VPA == "NO VPA exposure", 0, NA))
for (v in c("AgeOnset","AgeAtFirstScan","pBAG_winsorized"))
  meta[[v]] <- suppressWarnings(as.numeric(meta[[v]]))

cov <- read.delim(COV, stringsAsFactors = FALSE)
dat <- merge(meta[, c("FID","IID","pBAG_winsorized","AgeAtFirstScan","AgeOnset",
                      "Syndrome","TC_seizures","DRE_bin","VPA_bin")],
             cov, by = c("FID","IID"))

# 5 PRS at their step 4 best threshold, z-scored, |z| > 3 set to NA
m3 <- read.delim(M3, stringsAsFactors = FALSE)
best <- do.call(rbind, lapply(unique(m3$GWAS), function(g) {
  s <- m3[m3$GWAS == g, ]; s[which.min(s$P_value), c("GWAS","Threshold")]
}))
gwas_files <- list(
  "Epilepsy"        = "PRS_epilepsy_syndrome.all_score",
  "Intelligence"    = "PRS_intelligence_syndrome.all_score",
  "Longevity"       = "PRS_longevity_syndrome.all_score",
  "Frailty"         = "PRS_frailty_syndrome.all_score",
  "Brain Age (GWM)" = "PRS_brainage_gwm_syndrome.all_score"
)
prs_var <- c("Epilepsy"="PRS_epi","Intelligence"="PRS_intel","Longevity"="PRS_long",
             "Frailty"="PRS_frail","Brain Age (GWM)"="PRS_pBAG")

for (g in names(gwas_files)) {
  pt <- paste0("Pt_", best$Threshold[best$GWAS == g])
  sc <- read.table(file.path(SCO, gwas_files[[g]]), header = TRUE, stringsAsFactors = FALSE)
  if (!(pt %in% names(sc))) next
  tmp <- sc[, c("FID","IID", pt)]
  names(tmp)[3] <- prs_var[[g]]
  tmp[[prs_var[[g]]]] <- as.numeric(scale(tmp[[prs_var[[g]]]]))
  tmp[[prs_var[[g]]]][abs(tmp[[prs_var[[g]]]]) > 3] <- NA
  dat <- merge(dat, tmp, by = c("FID","IID"), all.x = TRUE)
}
cat("Merged n =", nrow(dat), "\n\n")

# Part A: tertile and quartile stratification
strata <- function(x, n) {
  q <- quantile(x, probs = seq(0, 1, length.out = n + 1), na.rm = TRUE)
  cut(x, breaks = q, include.lowest = TRUE, labels = paste0("Q", seq_len(n)))
}

strat_sum <- data.frame(); strat_plots <- list()
for (g in names(prs_var)) {
  v <- prs_var[[g]]
  df <- dat[, c("pBAG_winsorized", v)]
  df <- df[complete.cases(df), ]
  df$T3 <- strata(df[[v]], 3); df$Q4 <- strata(df[[v]], 4)

  aov_t <- summary(aov(pBAG_winsorized ~ T3, data = df))[[1]]
  aov_q <- summary(aov(pBAG_winsorized ~ Q4, data = df))[[1]]
  kw_t  <- kruskal.test(pBAG_winsorized ~ T3, data = df)
  kw_q  <- kruskal.test(pBAG_winsorized ~ Q4, data = df)
  lt_t  <- summary(lm(pBAG_winsorized ~ as.numeric(T3), data = df))$coefficients[2, ]
  lt_q  <- summary(lm(pBAG_winsorized ~ as.numeric(Q4), data = df))$coefficients[2, ]

  for (bin in list(list("Tertile", aov_t, kw_t, lt_t),
                   list("Quartile", aov_q, kw_q, lt_q))) {
    strat_sum <- rbind(strat_sum, data.frame(
      GWAS = g, Bin = bin[[1]], N = nrow(df),
      F_stat = round(bin[[2]][1, "F value"], 3),
      ANOVA_P = signif(bin[[2]][1, "Pr(>F)"], 3),
      KW_P = signif(bin[[3]]$p.value, 3),
      Trend_beta = round(bin[[4]]["Estimate"], 3),
      Trend_P = signif(bin[[4]]["Pr(>|t|)"], 3)
    ))
  }

  df_long <- rbind(
    data.frame(Bin = "Tertile", Stratum = df$T3, pBAG = df$pBAG_winsorized),
    data.frame(Bin = "Quartile", Stratum = df$Q4, pBAG = df$pBAG_winsorized))
  strat_plots[[g]] <- ggplot(df_long, aes(Stratum, pBAG, fill = Stratum)) +
    geom_boxplot(outlier.size = 0.7, alpha = 0.85) +
    geom_jitter(width = 0.15, size = 0.4, alpha = 0.3) +
    facet_wrap(~ Bin, scales = "free_x") +
    scale_fill_brewer(palette = "Blues") +
    labs(x = "", y = "pBAG (winsorized)",
         title = paste0("PRS - ", g),
         subtitle = sprintf("Tertile: ANOVA p=%.3g | trend p=%.3g  *  Quartile: ANOVA p=%.3g | trend p=%.3g",
                            aov_t[1,"Pr(>F)"], lt_t["Pr(>|t|)"],
                            aov_q[1,"Pr(>F)"], lt_q["Pr(>|t|)"])) +
    theme_minimal(base_size = 11) + theme(legend.position = "none")
}
cat("Part A summary:\n"); print(strat_sum, row.names = FALSE)
write.table(strat_sum, file.path(OUT, "PRS_stratification_summary.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Part B: PRS x clinical factor interactions
factors <- list(Syndrome = "Syndrome", DRE = "DRE_bin",
                Convulsions = "TC_seizures", VPA = "VPA_bin")
int_sum <- data.frame(); int_plots <- list()

for (g in names(prs_var)) {
  v <- prs_var[[g]]
  for (fn in names(factors)) {
    fv <- factors[[fn]]
    df <- dat[, c("pBAG_winsorized", v, fv, "PC1","PC2","PC3","PC4","Sex","AgeAtFirstScan")]
    df <- df[complete.cases(df), ]
    if (nrow(df) < 30) next
    df$FACTOR <- if (fn == "Syndrome") factor(df[[fv]], levels = c("Focal","Generalised"))
                 else factor(df[[fv]], levels = c(0, 1))
    df$PRS <- df[[v]]
    fit <- lm(pBAG_winsorized ~ PRS * FACTOR + PC1+PC2+PC3+PC4 + Sex + AgeAtFirstScan,
              data = df)
    co  <- summary(fit)$coefficients
    int_nm <- grep(":FACTOR", rownames(co), value = TRUE)
    if (length(int_nm) == 0) next
    ir <- co[int_nm, ]; mr <- co["PRS", ]
    lvs <- levels(df$FACTOR)
    slopes <- lapply(lvs, function(lv) {
      d <- df[df$FACTOR == lv, ]
      f <- lm(pBAG_winsorized ~ PRS + PC1+PC2+PC3+PC4 + Sex + AgeAtFirstScan, data = d)
      summary(f)$coefficients["PRS", ]
    })

    int_sum <- rbind(int_sum, data.frame(
      GWAS = g, Factor = fn, N_total = nrow(df),
      N_lvl0 = sum(df$FACTOR == lvs[1]), N_lvl1 = sum(df$FACTOR == lvs[2]),
      PRS_main_beta = round(mr["Estimate"], 3),
      PRS_main_P = signif(mr["Pr(>|t|)"], 3),
      Interaction_beta = round(ir["Estimate"], 3),
      Interaction_SE = round(ir["Std. Error"], 3),
      Interaction_P = signif(ir["Pr(>|t|)"], 3),
      Slope_lvl0 = round(slopes[[1]]["Estimate"], 3),
      P_lvl0 = signif(slopes[[1]]["Pr(>|t|)"], 3),
      Slope_lvl1 = round(slopes[[2]]["Estimate"], 3),
      P_lvl1 = signif(slopes[[2]]["Pr(>|t|)"], 3)
    ))

    df$FACTOR_label <- factor(df$FACTOR, levels = lvs,
      labels = if (fn == "Syndrome") lvs else c(paste0("No ", fn), fn))
    int_plots[[paste(g, fn, sep = "_")]] <- ggplot(df,
        aes(PRS, pBAG_winsorized, color = FACTOR_label)) +
      geom_point(alpha = 0.25, size = 0.9) +
      geom_smooth(method = "lm", se = TRUE, linewidth = 1) +
      scale_color_manual(values = c("#2166AC","#B2182B"), name = "") +
      labs(x = paste0(g, " PRS (z)"), y = "pBAG (winsorized)",
           title = paste0(g, " x ", fn),
           subtitle = sprintf("Int beta=%.2f p=%.3g; slope %s=%.2f (p=%.3g); slope %s=%.2f (p=%.3g)",
              ir["Estimate"], ir["Pr(>|t|)"],
              lvs[1], slopes[[1]]["Estimate"], slopes[[1]]["Pr(>|t|)"],
              lvs[2], slopes[[2]]["Estimate"], slopes[[2]]["Pr(>|t|)"])) +
      theme_minimal(base_size = 10) +
      theme(legend.position = "bottom",
            plot.subtitle = element_text(size = 7))
  }
}
cat("\nPart B summary:\n"); print(int_sum, row.names = FALSE)
write.table(int_sum, file.path(OUT, "PRS_interaction_summary.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Heatmap of interaction p-values (5 PRS x 4 factors)
ih <- int_sum[, c("GWAS","Factor","Interaction_beta","Interaction_P")]
ih$sig <- ifelse(ih$Interaction_P < 0.001, "***",
          ifelse(ih$Interaction_P < 0.01, "**",
          ifelse(ih$Interaction_P < 0.05, "*",
          ifelse(ih$Interaction_P < 0.1, ".", ""))))
ih$GWAS   <- factor(ih$GWAS, levels = c("Epilepsy","Intelligence","Longevity","Frailty","Brain Age (GWM)"))
ih$Factor <- factor(ih$Factor, levels = c("Syndrome","DRE","Convulsions","VPA"))

p_heat <- ggplot(ih, aes(Factor, GWAS, fill = Interaction_beta)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = sprintf("beta=%.2f\np=%.3g%s",
                                Interaction_beta, Interaction_P, sig)), size = 3) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, name = "Interaction\nbeta") +
  labs(title = "PRS x clinical factor interaction on pBAG",
       x = "", y = "") +
  theme_minimal(base_size = 12) + theme(panel.grid = element_blank())

pdf(file.path(OUT, "PRS_stratification.pdf"), width = 12, height = 8)
print(ggtexttable(strat_sum, rows = NULL, theme = ttheme("light")))
for (g in names(strat_plots)) print(strat_plots[[g]])
print(p_heat)
for (g in names(prs_var)) {
  keep <- startsWith(names(int_plots), paste0(g, "_"))
  ps <- int_plots[keep]
  if (length(ps) == 0) next
  print(ggarrange(plotlist = ps, ncol = 2, nrow = 2,
                  common.legend = TRUE, legend = "bottom"))
}
dev.off()
cat("Saved PDF:", file.path(OUT, "PRS_stratification.pdf"), "\n")
