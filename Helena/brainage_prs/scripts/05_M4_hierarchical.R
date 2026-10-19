#!/usr/bin/env Rscript
# Hierarchical block regression of pBAG: PCs, + demographics, + epilepsy factors, + the 5 PRS.
# Usage: Rscript scripts/05_M4_hierarchical.R config.sh (needs PRS_pBAG_linear_results.txt from step 4)

suppressPackageStartupMessages({
  library(ggplot2); library(ggpubr); library(reshape2)
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

names(meta)[names(meta) == "Age at onset"]            <- "AgeOnset"
names(meta)[names(meta) == "Years.Of.Epilepsy"]       <- "Duration"
names(meta)[names(meta) == "EpilepsySyndrome_forPRS"] <- "Syndrome"
names(meta)[names(meta) == "HistoryOfConvulsions"]    <- "TC_seizures"

meta$Syndrome    <- ifelse(meta$Syndrome %in% c("Focal","Generalised"), meta$Syndrome, NA)
meta$TC_seizures <- ifelse(meta$TC_seizures == "History of convulsions", 1,
                    ifelse(meta$TC_seizures == "NO History of convulsions", 0, NA))
meta$VPA_bin     <- ifelse(meta$VPA == "VPA exposure", 1,
                    ifelse(meta$VPA == "NO VPA exposure", 0, NA))
for (v in c("AgeOnset","Duration","AgeAtFirstScan","pBAG_winsorized"))
  meta[[v]] <- suppressWarnings(as.numeric(meta[[v]]))

cov <- read.delim(COV, stringsAsFactors = FALSE)
dat <- merge(meta[, c("FID","IID","pBAG_winsorized","AgeAtFirstScan",
                      "AgeOnset","Duration","Syndrome","TC_seizures","VPA_bin")],
             cov, by = c("FID","IID"))

m3 <- read.delim(M3, stringsAsFactors = FALSE)
best <- do.call(rbind, lapply(unique(m3$GWAS), function(g) {
  s <- m3[m3$GWAS == g, ]; s[which.min(s$P_value), c("GWAS","Threshold")]
}))
cat("Best Pt per GWAS:\n"); print(best, row.names = FALSE)

gwas_files <- list(
  "Epilepsy"        = "PRS_epilepsy_syndrome.all_score",
  "Intelligence"    = "PRS_intelligence_syndrome.all_score",
  "Longevity"       = "PRS_longevity_syndrome.all_score",
  "Frailty"         = "PRS_frailty_syndrome.all_score",
  "Brain Age (GWM)" = "PRS_brainage_gwm_syndrome.all_score"
)
prs_var <- c("Epilepsy"="PRS_epi", "Intelligence"="PRS_intel",
             "Longevity"="PRS_long", "Frailty"="PRS_frail",
             "Brain Age (GWM)"="PRS_pBAG")

for (g in names(gwas_files)) {
  pt <- paste0("Pt_", best$Threshold[best$GWAS == g])
  sc <- read.table(file.path(SCO, gwas_files[[g]]), header = TRUE, stringsAsFactors = FALSE)
  if (!(pt %in% names(sc))) next
  tmp <- sc[, c("FID","IID", pt)]
  names(tmp)[3] <- prs_var[[g]]
  tmp[[prs_var[[g]]]] <- as.numeric(scale(tmp[[prs_var[[g]]]]))     # per-SD
  tmp[[prs_var[[g]]]][abs(tmp[[prs_var[[g]]]]) > 3] <- NA            # clip +-3 SD
  dat <- merge(dat, tmp, by = c("FID","IID"), all.x = TRUE)
}

# Duration (Years.Of.Epilepsy) is AgeAtScan - AgeOnset, perfectly collinear with
# blocks 1-2, so it is read but kept out of the models.
mvars <- c("pBAG_winsorized","AgeAtFirstScan","Sex","AgeOnset","Syndrome",
           "TC_seizures","VPA_bin","PC1","PC2","PC3","PC4", as.character(prs_var))
dat_cc <- dat[complete.cases(dat[, mvars]), ]
cat("\nComplete-case n =", nrow(dat_cc), "\n\n")

f0 <- pBAG_winsorized ~ PC1+PC2+PC3+PC4
f1 <- update(f0, . ~ . + AgeAtFirstScan + Sex)
f2 <- update(f1, . ~ . + AgeOnset + Syndrome + TC_seizures + VPA_bin)
f3 <- update(f2, as.formula(paste(". ~ . +", paste(prs_var, collapse = " + "))))
m0 <- lm(f0, data = dat_cc); m1 <- lm(f1, data = dat_cc)
m2 <- lm(f2, data = dat_cc); m3lm <- lm(f3, data = dat_cc)

r2 <- c("Covariates only (PCs)" = summary(m0)$r.squared,
        "Block 1: + Demographics" = summary(m1)$r.squared,
        "Block 2: + Epilepsy factors" = summary(m2)$r.squared,
        "Block 3: + PRS (genetics)" = summary(m3lm)$r.squared)
an <- anova(m0, m1, m2, m3lm); cat("ANOVA:\n"); print(an)

blocks <- data.frame(
  Block = names(r2),
  R2 = round(as.numeric(r2), 5),
  Adj_R2 = round(c(summary(m0)$adj.r.squared, summary(m1)$adj.r.squared,
                   summary(m2)$adj.r.squared, summary(m3lm)$adj.r.squared), 5),
  Delta_R2 = round(c(NA, diff(as.numeric(r2))), 5),
  F  = round(c(NA, an$F[2:4]), 3),
  df1 = c(NA, an$Df[2:4]),
  df2 = c(NA, an$Res.Df[2:4]),
  P  = format(c(NA, an$`Pr(>F)`[2:4]), digits = 3, scientific = TRUE)
)
write.table(blocks, file.path(OUT, "hierarchical_regression_blocks.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("\nBlock-wise:\n"); print(blocks, row.names = FALSE)

co <- as.data.frame(summary(m3lm)$coefficients)
co$Predictor <- rownames(co); co <- co[, c(5,1,2,3,4)]
names(co) <- c("Predictor","Estimate","SE","t","P")
write.table(co, file.path(OUT, "hierarchical_regression_fullmodel.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

plot_df <- data.frame(Block = factor(names(r2), levels = names(r2)),
                      R2 = as.numeric(r2),
                      Delta_R2 = c(as.numeric(r2)[1], diff(as.numeric(r2))))
plot_long <- reshape2::melt(plot_df, id.vars = "Block",
                            variable.name = "Type", value.name = "R2")
plot_long$Type <- factor(plot_long$Type, levels = c("R2","Delta_R2"),
                         labels = c("Cumulative R2","Incremental dR2"))

pA <- ggplot(plot_long, aes(Block, R2, fill = Type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = sprintf("%.3f", R2)),
            position = position_dodge(width = 0.8), vjust = -0.4, size = 3.5) +
  scale_fill_manual(values = c("Cumulative R2"="#2166AC","Incremental dR2"="#B2182B"),
                    name = "") +
  labs(x = "", y = expression(R^2),
       title = "Hierarchical regression of pBAG (winsorized)",
       subtitle = sprintf("n = %d (complete-case)", nrow(dat_cc))) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 18, hjust = 1), legend.position = "top")

tab <- blocks
tab$Delta_R2 <- ifelse(is.na(tab$Delta_R2), "-", sprintf("%.4f", tab$Delta_R2))
tab$F        <- ifelse(is.na(tab$F), "-", sprintf("%.2f", tab$F))
tab$P        <- ifelse(grepl("NA", tab$P), "-", tab$P)
pB <- ggtexttable(tab[, c("Block","R2","Adj_R2","Delta_R2","F","df1","df2","P")],
                  rows = NULL, theme = ttheme("light"))

cp <- co[!(co$Predictor %in% c("(Intercept)","PC1","PC2","PC3","PC4")), ]
lbl <- c(AgeAtFirstScan = "Age at scan", Sex = "Sex",
         AgeOnset = "Age at onset", Duration = "Duration",
         SyndromeGeneralised = "Syndrome (Gen vs Foc)",
         TC_seizures = "History of convulsions", VPA_bin = "VPA exposure",
         PRS_epi = "PRS - Epilepsy", PRS_intel = "PRS - Intelligence",
         PRS_long = "PRS - Longevity", PRS_frail = "PRS - Frailty",
         PRS_pBAG = "PRS - Brain Age (GWM)")
cp$Label <- ifelse(cp$Predictor %in% names(lbl), lbl[cp$Predictor], cp$Predictor)
cp$ci_lo <- cp$Estimate - 1.96*cp$SE
cp$ci_hi <- cp$Estimate + 1.96*cp$SE
cp$Block <- ifelse(cp$Predictor %in% c("AgeAtFirstScan","Sex"), "1. Demographics",
            ifelse(cp$Predictor %in% c("AgeOnset","Duration","SyndromeGeneralised",
                                       "TC_seizures","VPA_bin"),
                   "2. Epilepsy factors", "3. Genetics (PRS)"))
cp$Label <- factor(cp$Label,
                   levels = cp$Label[order(cp$Block, cp$Estimate)])
cp$Sig <- factor(ifelse(cp$P < 0.05, "p < 0.05",
                 ifelse(cp$P < 0.10, "p < 0.1", "n.s.")),
                 levels = c("p < 0.05","p < 0.1","n.s."))

pC <- ggplot(cp, aes(Estimate, Label, color = Sig)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.25, linewidth = 0.7) +
  geom_point(size = 2.8) +
  facet_grid(Block ~ ., scales = "free_y", space = "free_y") +
  scale_color_manual(values = c("p < 0.05"="#B2182B","p < 0.1"="#F4A582","n.s."="#636363"),
                     name = "") +
  labs(title = "Block 3 coefficients (Demographics + Epilepsy + PRS)",
       subtitle = "95% CI; PRS are z-scored (per-SD effect)",
       x = "Coefficient (pBAG units)", y = "") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom",
        strip.text.y = element_text(angle = 0, face = "bold"))

pdf(file.path(OUT, "hierarchical_regression_pBAG.pdf"), width = 11, height = 8)
print(pA); print(ggarrange(pB, ncol = 1)); print(pC)
dev.off()
cat("Saved PDF:", file.path(OUT, "hierarchical_regression_pBAG.pdf"), "\n")
