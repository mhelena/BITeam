#!/usr/bin/env Rscript
# Clinical-only baseline for pBAG (no PRS): PCs, + age at scan and sex, + onset, syndrome, TC seizures, VPA.
# Usage: Rscript scripts/03_M2_clinical_only.R config.sh

suppressPackageStartupMessages(library(ggplot2))

args <- commandArgs(trailingOnly = FALSE)
source(file.path(dirname(sub("^--file=", "", grep("^--file=", args, value = TRUE))), "config.R"))
cfg <- read_config()

META <- cfg("METADATA")
COV  <- cfg("COVARIATES")
OUT  <- cfg("OUT_DIR")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

meta <- read.delim(META, stringsAsFactors = FALSE, check.names = FALSE)
meta <- meta[!(meta$MRI_findings %in% c("Post surgery", "Large lesion")), ]

names(meta)[names(meta) == "Age at onset"]            <- "AgeOnset"
names(meta)[names(meta) == "EpilepsySyndrome_forPRS"] <- "Syndrome"
names(meta)[names(meta) == "HistoryOfConvulsions"]    <- "TC_seizures"

meta$Syndrome    <- ifelse(meta$Syndrome %in% c("Focal","Generalised"), meta$Syndrome, NA)
meta$TC_seizures <- ifelse(meta$TC_seizures == "History of convulsions", 1,
                    ifelse(meta$TC_seizures == "NO History of convulsions", 0, NA))
meta$VPA_bin     <- ifelse(meta$VPA == "VPA exposure", 1,
                    ifelse(meta$VPA == "NO VPA exposure", 0, NA))
for (v in c("AgeOnset","AgeAtFirstScan","pBAG_winsorized"))
  meta[[v]] <- suppressWarnings(as.numeric(meta[[v]]))

cov <- read.delim(COV, stringsAsFactors = FALSE)
dat <- merge(meta[, c("FID","IID","pBAG_winsorized","AgeAtFirstScan",
                      "AgeOnset","Syndrome","TC_seizures","VPA_bin")],
             cov, by = c("FID","IID"))
dat <- dat[complete.cases(dat[, c("pBAG_winsorized","AgeAtFirstScan","Sex",
                                   "AgeOnset","Syndrome","TC_seizures","VPA_bin",
                                   paste0("PC", 1:4))]), ]
cat("n =", nrow(dat), "\n")

m0 <- lm(pBAG_winsorized ~ PC1+PC2+PC3+PC4, data = dat)
m1 <- update(m0, . ~ . + AgeAtFirstScan + Sex)
m2 <- update(m1, . ~ . + AgeOnset + Syndrome + TC_seizures + VPA_bin)
an <- anova(m0, m1, m2); print(an)

cat(sprintf("\nR2: PCs %.4f | +demo %.4f (dR2 %.4f, p=%.3g) | +clin %.4f (dR2 %.4f, p=%.3g)\n",
    summary(m0)$r.squared,
    summary(m1)$r.squared, summary(m1)$r.squared - summary(m0)$r.squared, an$`Pr(>F)`[2],
    summary(m2)$r.squared, summary(m2)$r.squared - summary(m1)$r.squared, an$`Pr(>F)`[3]))

co <- as.data.frame(summary(m2)$coefficients)
co$Predictor <- rownames(co); names(co) <- c("Estimate","SE","t","P","Predictor")
co <- co[, c("Predictor","Estimate","SE","t","P")]
write.table(co, file.path(OUT, "clinical_only_coefs.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

pl <- co[!(co$Predictor %in% c("(Intercept)","PC1","PC2","PC3","PC4")), ]
pl$Label <- c(AgeAtFirstScan = "Age at scan", Sex = "Sex (M vs F)",
              AgeOnset = "Age at onset",
              SyndromeGeneralised = "Generalised vs Focal",
              TC_seizures = "History of TC seizures",
              VPA_bin = "Current VPA use")[pl$Predictor]
pl$ci_lo <- pl$Estimate - 1.96*pl$SE
pl$ci_hi <- pl$Estimate + 1.96*pl$SE
pl$sig   <- ifelse(pl$P < 0.05, "p < 0.05",
            ifelse(pl$P < 0.10, "marginal", "ns"))
pl <- pl[order(-abs(pl$Estimate)), ]
pl$Label <- factor(pl$Label, levels = rev(pl$Label))

p <- ggplot(pl, aes(Estimate, Label, color = sig)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.25, linewidth = 0.7) +
  geom_point(size = 3.2) +
  scale_color_manual(values = c("p < 0.05"="#C0392B","marginal"="#B86E10","ns"="#1F3A68")) +
  labs(x = "Effect on pBAG (years)", y = NULL, color = NULL,
       title = "Clinical and demographic predictors of pBAG",
       subtitle = sprintf("n = %d  |  R\u00b2 = %.3f",
                          nrow(dat), summary(m2)$r.squared)) +
  theme_minimal(base_size = 13) +
  theme(panel.grid.major.y = element_blank(),
        plot.title = element_text(face = "bold", color = "#1F3A68"),
        legend.position = "bottom")
ggsave(file.path(OUT, "clinical_only_forest.png"), p,
       width = 9, height = 4.5, dpi = 200)
cat("Saved: clinical_only_coefs.txt + clinical_only_forest.png\n")
