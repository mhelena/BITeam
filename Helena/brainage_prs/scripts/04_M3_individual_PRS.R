#!/usr/bin/env Rscript
# Regresses pBAG on each of the 5 selected PRS at all 7 PRSice thresholds (pBAG ~ PRS + PC1-4 + Sex).
# Usage: Rscript scripts/04_M3_individual_PRS.R config.sh

library(ggplot2)
library(ggpubr)

args <- commandArgs(trailingOnly = FALSE)
source(file.path(dirname(sub("^--file=", "", grep("^--file=", args, value = TRUE))), "config.R"))
cfg <- read_config()

META <- cfg("METADATA")
COV  <- cfg("COVARIATES")
SCO  <- cfg("PRS_DIR")
OUT  <- cfg("OUT_DIR")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

meta <- read.delim(META, stringsAsFactors = FALSE, check.names = FALSE)
meta <- meta[!(meta$MRI_findings %in% c("Post surgery", "Large lesion")), ]
cat("After exclusion:", nrow(meta), "\n")

cov <- read.delim(COV, stringsAsFactors = FALSE)
dat <- merge(meta[, c("FID","IID","pBAG_winsorized")], cov, by = c("FID","IID"))
dat$pBAG_winsorized <- suppressWarnings(as.numeric(dat$pBAG_winsorized))

# The 5 PRS kept after LDSC genetic correlation screening of the 8 bases
gwas <- list(
  epilepsy     = list(label="Epilepsy",        file="PRS_epilepsy_syndrome.all_score"),
  intelligence = list(label="Intelligence",    file="PRS_intelligence_syndrome.all_score"),
  longevity    = list(label="Longevity",       file="PRS_longevity_syndrome.all_score"),
  frailty      = list(label="Frailty",         file="PRS_frailty_syndrome.all_score"),
  brainage_gwm = list(label="Brain Age (GWM)", file="PRS_brainage_gwm_syndrome.all_score")
)
thr <- c("Pt_0.0001","Pt_0.001","Pt_0.01","Pt_0.05","Pt_0.1","Pt_0.5","Pt_1")
thr_lab <- c("1e-04","0.001","0.01","0.05","0.1","0.5","1")

res <- data.frame()
for (g in names(gwas)) {
  info <- gwas[[g]]
  sc <- read.table(file.path(SCO, info$file), header = TRUE, stringsAsFactors = FALSE)
  df <- merge(dat, sc, by = c("FID","IID"))
  null <- lm(pBAG_winsorized ~ PC1+PC2+PC3+PC4+Sex, data = df)
  null_r2 <- summary(null)$r.squared

  for (i in seq_along(thr)) {
    df$PRS <- df[[thr[i]]]
    # drop +-3 SD outliers
    m <- mean(df$PRS); s <- sd(df$PRS)
    d <- df[abs(df$PRS - m) <= 3*s, ]
    full <- lm(pBAG_winsorized ~ PRS + PC1+PC2+PC3+PC4+Sex, data = d)
    co <- summary(full)$coefficients["PRS", ]
    res <- rbind(res, data.frame(
      GWAS = info$label, Threshold = thr_lab[i],
      N = nrow(d), Outliers_removed = nrow(df) - nrow(d),
      Null_R2 = round(null_r2, 6),
      Full_R2 = round(summary(full)$r.squared, 6),
      PRS_R2  = round(summary(full)$r.squared - null_r2, 6),
      Coefficient = co["Estimate"], SE = co["Std. Error"],
      T_value = co["t value"], P_value = co["Pr(>|t|)"]
    ))
  }
}
res <- res[order(res$P_value), ]; rownames(res) <- NULL
outfile <- file.path(OUT, "PRS_pBAG_linear_results.txt")
write.table(res, outfile, sep = "\t", row.names = FALSE, quote = FALSE)
cat("Saved:", outfile, "\n\n")

cat("Best threshold per GWAS:\n")
best <- do.call(rbind, lapply(unique(res$GWAS), function(g) {
  s <- res[res$GWAS == g, ]; s[which.min(s$P_value), ]
}))
for (i in seq_len(nrow(best))) {
  b <- best[i, ]
  cat(sprintf("  %-18s Pt=%s  R2=%.4f  beta=%.3f  P=%.2e%s\n",
              b$GWAS, b$Threshold, b$PRS_R2, b$Coefficient, b$P_value,
              ifelse(b$P_value < 0.05, " *", "")))
}

best$Sig <- ifelse(best$P_value < 0.05, "p < 0.05", "n.s.")
best$GWAS <- factor(best$GWAS, levels = best$GWAS[order(best$PRS_R2)])

p_bar <- ggplot(best, aes(GWAS, PRS_R2, fill = Sig)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_fill_manual(values = c("#0072B2", "#999999"), name = "") +
  coord_flip() +
  labs(x = "", y = expression(paste("PRS ", R^2, " (incremental)")),
       title = "PRS association with pBAG (winsorized)",
       subtitle = "Best p-value threshold per GWAS | adjusted for PCs + Sex") +
  theme_minimal(base_size = 14)

scatters <- list()
for (g in names(gwas)) {
  info <- gwas[[g]]
  sc <- read.table(file.path(SCO, info$file), header = TRUE, stringsAsFactors = FALSE)
  df <- merge(dat, sc, by = c("FID","IID"))
  b <- best[best$GWAS == info$label, ]
  df$PRS <- df[[paste0("Pt_", b$Threshold)]]
  m <- mean(df$PRS); s <- sd(df$PRS); df <- df[abs(df$PRS - m) <= 3*s, ]
  df$PRS_scaled <- scale(df$PRS)[, 1]
  scatters[[length(scatters)+1]] <- ggplot(df, aes(PRS_scaled, pBAG_winsorized)) +
    geom_point(alpha = 0.3, size = 1.5, color = "#0072B2") +
    geom_smooth(method = "lm", color = "red", se = TRUE) +
    labs(x = paste0("PRS (z) - Pt=", b$Threshold), y = "pBAG (winsorized)",
         title = info$label,
         subtitle = sprintf("R\u00b2=%.4f, P=%.2e", b$PRS_R2, b$P_value)) +
    theme_minimal(base_size = 12)
}

pdf(file.path(OUT, "PRS_pBAG_linear.pdf"), width = 14, height = 10)
print(p_bar)
print(ggarrange(plotlist = scatters, ncol = 3, nrow = 2))
dev.off()
cat("Saved PDF:", file.path(OUT, "PRS_pBAG_linear.pdf"), "\n")
