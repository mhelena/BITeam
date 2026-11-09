# Builds the simulated case-only cohort used by the notebooks (Source it once, about 10 seconds).
# Writes data/toy/: genotypes, clinical table, genetic findings, two base GWAS and the simulated truth.
# Nobody real is in these data.

source("config.R")
set.seed(2026)
dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)

N_POP    <- 40000    # population simulated
BLOCK    <- 10       # SNPs per LD block
N_BLOCKS <- 1500     # 15,000 SNPs in total
N_CAUSAL <- 200      # blocks carrying a causal epilepsy variant
H2       <- 0.4      # liability variance explained by common variants
PREV     <- 0.03     # proportion of the population who develop epilepsy
F_ST     <- 0.02     # allele frequency difference between two ancestry groups
PROP_B   <- 0.15     # share of ancestry group B
N_GWAS   <- c(epilepsy = 5000, height = 100000)   # base GWAS sample sizes

# Rare variants: population frequency and effect on liability (SD units)
MONO <- data.frame(
  category  = c("high_PTV", "high_missense", "incomplete", "VUS"),
  freq_A    = c(0.002,      0.002,           0.008,        0.04),
  effect    = c(3.3,        2.3,             1.5,          0)
)
FOUNDER_NPRL3_B <- 0.02     # extra NPRL3 carriers in ancestry group B (a founder variant)

# Genome layout
chr_len <- c(248, 242, 198, 190, 181, 171, 159, 145, 138, 134, 135, 133, 114, 107, 102, 90, 83, 80, 59, 64, 47, 51) * 1e6
blocks_per_chr <- pmax(1, round(N_BLOCKS * chr_len / sum(chr_len)))
blocks_per_chr[1] <- blocks_per_chr[1] + N_BLOCKS - sum(blocks_per_chr)
blk <- data.frame(block = seq_len(N_BLOCKS), chr = rep(1:22, blocks_per_chr))
blk$start <- unlist(lapply(1:22, function(c)
  sort(sample(seq(1e6, chr_len[c] - 1e6, by = 4e5), blocks_per_chr[c]))))
blk$rho <- runif(N_BLOCKS, 0.3, 0.85)                 # LD strength inside each block

snp <- data.frame(block = rep(blk$block, each = BLOCK), j = rep(1:BLOCK, N_BLOCKS))
snp$chr <- blk$chr[snp$block]
snp$bp  <- blk$start[snp$block] + (snp$j - 1) * 5000 + sample(0:4000, nrow(snp), TRUE)
snp$id  <- paste0("rs", sample(1e6:9e7, nrow(snp)))
M <- nrow(snp)

# alleles: A1 is the allele counted in the genotypes; ~3% ambiguous (A/T, C/G)
pairs_ok  <- list(c("A", "C"), c("A", "G"), c("C", "T"), c("G", "T"), c("C", "A"), c("T", "G"))
pairs_amb <- list(c("A", "T"), c("C", "G"))
al <- t(sapply(seq_len(M), function(i) if (runif(1) < 0.03) sample(pairs_amb, 1)[[1]] else sample(pairs_ok, 1)[[1]]))
snp$A1 <- al[, 1]; snp$A2 <- al[, 2]

# allele frequencies in the two ancestry groups (Balding-Nichols)
p_anc <- runif(M, 0.03, 0.5)
bn <- function(p) pmin(pmax(rbeta(length(p), p * (1 - F_ST) / F_ST, (1 - p) * (1 - F_ST) / F_ST), 0.005), 0.995)
p_pop <- rbind(A = bn(p_anc), B = bn(p_anc))          # 2 x M

# simulate one LD block for a set of people (pop = 1 for A, 2 for B)
sim_block <- function(b, pop) {
  cols <- which(snp$block == b)
  n <- length(pop)
  thr <- qnorm(p_pop[pop, cols, drop = FALSE])
  g <- matrix(0L, n, length(cols))
  for (h in 1:2) {
    z <- sqrt(blk$rho[b]) * rnorm(n) + sqrt(1 - blk$rho[b]) * matrix(rnorm(n * length(cols)), n)
    g <- g + (z < thr)
  }
  g
}

# 1. Population and liability
message("Simulating the population ...")
pop <- ifelse(runif(N_POP) < PROP_B, 2L, 1L)
causal_blocks <- sort(sample(N_BLOCKS, N_CAUSAL))
height_blocks <- sort(sample(setdiff(seq_len(N_BLOCKS), causal_blocks), N_CAUSAL))
causal_snp <- sapply(causal_blocks, function(b) which(snp$block == b)[sample(BLOCK, 1)])
height_snp <- sapply(height_blocks, function(b) which(snp$block == b)[sample(BLOCK, 1)])
p_pool <- (1 - PROP_B) * p_pop["A", ] + PROP_B * p_pop["B", ]
beta_std <- rnorm(N_CAUSAL, 0, sqrt(H2 / N_CAUSAL))

pop_geno <- vector("list", N_CAUSAL)                  # keep causal blocks for everyone
G <- numeric(N_POP)
for (k in seq_along(causal_blocks)) {
  gb <- sim_block(causal_blocks[k], pop)
  pop_geno[[k]] <- gb
  s <- causal_snp[k]
  jj <- which(which(snp$block == causal_blocks[k]) == s)
  G <- G + beta_std[k] * (gb[, jj] - 2 * p_pool[s]) / sqrt(2 * p_pool[s] * (1 - p_pool[s]))
}
scale_g <- sqrt(H2 / var(G)); G <- G * scale_g; beta_std <- beta_std * scale_g

# rare variants
cat_prob <- rbind(A = c(MONO$freq_A, 1 - sum(MONO$freq_A)),
                  B = c(MONO$freq_A[1:2], MONO$freq_A[3] + FOUNDER_NPRL3_B, MONO$freq_A[4],
                        1 - sum(MONO$freq_A) - FOUNDER_NPRL3_B))
mono_cat <- character(N_POP)
for (p in 1:2) {
  idx <- which(pop == p)
  mono_cat[idx] <- sample(c(MONO$category, "none"), length(idx), TRUE, prob = cat_prob[p, ])
}
mono_eff <- c(setNames(MONO$effect, MONO$category), none = 0)[mono_cat]
L <- G + mono_eff + rnorm(N_POP, 0, sqrt(1 - H2))
case <- which(L > quantile(L, 1 - PREV))
n <- length(case)
message(sprintf("Cohort: %d people with epilepsy", n))

# 2. Cohort genotypes
message("Simulating cohort genotypes ...")
geno <- matrix(0L, n, M)
for (b in seq_len(N_BLOCKS)) {
  cols <- which(snp$block == b)
  k <- match(b, causal_blocks)
  geno[, cols] <- if (!is.na(k)) pop_geno[[k]][case, ] else sim_block(b, pop[case])
}
rm(pop_geno); invisible(gc())
geno_clean <- geno                                    # used later for the GWAS LD

iid <- sprintf("EPI%04d", sample(1000:9999, n))
cpop <- pop[case]; cmono <- mono_cat[case]

# 3. Problems for the QC session
qc_truth <- data.frame(IID = character(), problem = character())
# a duplicated sample (same person sequenced twice)
dup_from <- sample(n, 1)
geno <- rbind(geno, geno[dup_from, ])
iid <- c(iid, paste0(iid[dup_from], "_R"))
cpop <- c(cpop, cpop[dup_from]); cmono <- c(cmono, cmono[dup_from])
qc_truth <- rbind(qc_truth, data.frame(IID = iid[length(iid)], problem = paste("duplicate of", iid[dup_from])))
# three pairs of first-degree relatives (share half their genome)
rel <- matrix(sample(seq_len(n)[-dup_from], 6), ncol = 2)
for (r in seq_len(nrow(rel))) {
  shared <- which(snp$block %in% sample(N_BLOCKS, N_BLOCKS / 2))
  geno[rel[r, 2], shared] <- geno[rel[r, 1], shared]
  qc_truth <- rbind(qc_truth, data.frame(IID = iid[rel[r, 2]], problem = paste("relative of", iid[rel[r, 1]])))
}
n <- nrow(geno)
# samples with many missing genotypes
bad_miss <- sample(setdiff(seq_len(n), c(rel, dup_from, n)), 6)
for (i in bad_miss) geno[i, sample(M, round(0.12 * M))] <- NA
qc_truth <- rbind(qc_truth, data.frame(IID = iid[bad_miss], problem = "high missingness"))
# contaminated samples (too many heterozygous calls)
bad_het <- sample(setdiff(seq_len(n), c(rel, dup_from, n, bad_miss)), 3)
for (i in bad_het) { hom <- which(geno[i, ] %in% c(0L, 2L)); geno[i, sample(hom, round(0.3 * length(hom)))] <- 1L }
qc_truth <- rbind(qc_truth, data.frame(IID = iid[bad_het], problem = "excess heterozygosity"))
# poorly genotyped SNPs and SNPs failing Hardy-Weinberg
bad_snp_miss <- sample(M, 150)
for (s in bad_snp_miss) geno[sample(n, round(0.08 * n)), s] <- NA
bad_hwe <- sample(setdiff(seq_len(M), bad_snp_miss), 40)
for (s in bad_hwe) geno[sample(n, round(0.6 * n)), s] <- 1L

# 4. Write PLINK files
message("Writing PLINK files ...")
write_bed <- function(g, prefix) {
  # PLINK 1 .bed, SNP-major: 00 = two copies of A1, 10 = one, 11 = none, 01 = missing
  code <- ifelse(is.na(g), 1L, c(3L, 2L, 0L)[g + 1L])
  n <- nrow(code); pad <- (4 - n %% 4) %% 4
  if (pad) code <- rbind(code, matrix(0L, pad, ncol(code)))
  arr <- array(code, c(4, nrow(code) / 4, ncol(code)))
  bytes <- arr[1, , ] + 4L * arr[2, , ] + 16L * arr[3, , ] + 64L * arr[4, , ]
  con <- file(paste0(prefix, ".bed"), "wb")
  writeBin(as.raw(c(0x6c, 0x1b, 0x01)), con); writeBin(as.raw(as.vector(bytes)), con)
  close(con)
}
write_bed(geno, file.path(DATA_DIR, "toy_cohort"))
write.table(data.frame(snp$chr, snp$id, 0, snp$bp, snp$A1, snp$A2), file.path(DATA_DIR, "toy_cohort.bim"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(data.frame(iid, iid, 0, 0, 0, -9), file.path(DATA_DIR, "toy_cohort.fam"),
            sep = " ", quote = FALSE, row.names = FALSE, col.names = FALSE)

# 5. Clinical table and genetic findings
message("Writing clinical tables ...")
genes_high <- GENES_HIGH_PENETRANCE; genes_inc <- GENES_INCOMPLETE
gene <- ifelse(cmono %in% c("high_PTV", "high_missense"), sample(genes_high, n, TRUE),
        ifelse(cmono == "incomplete", ifelse(cpop == 2, "NPRL3", sample(genes_inc, n, TRUE)),
        ifelse(cmono == "VUS", sample(c(genes_high, genes_inc, "GRIN2A", "GABRB3", "CHD2"), n, TRUE), NA)))
consequence <- ifelse(cmono == "high_PTV", sample(c("stop_gained", "frameshift_variant", "splice_donor_variant"), n, TRUE),
               ifelse(cmono == "high_missense", "missense_variant",
               ifelse(cmono == "incomplete", sample(c("stop_gained", "frameshift_variant", "missense_variant"), n, TRUE,
                                                    prob = c(0.4, 0.4, 0.2)),
               ifelse(cmono == "VUS", "missense_variant", NA))))
acmg <- ifelse(cmono == "VUS", "VUS", ifelse(is.na(gene), NA, sample(c("Pathogenic", "Likely pathogenic"), n, TRUE)))
findings <- data.frame(IID = iid, gene = gene, variant = paste0(gene, "_toyvar", sample(10:99, n, TRUE)),
                       consequence = consequence, acmg = acmg)[!is.na(gene), ]
# a few people also carry an extra VUS (a second row for the same person)
extra <- sample(findings$IID[findings$acmg != "VUS"], 8)
findings <- rbind(findings, data.frame(IID = extra, gene = sample(c("GRIN2A", "CHD2", "GABRB3"), 8, TRUE),
                                       variant = paste0("extra_toyvar", 1:8), consequence = "missense_variant", acmg = "VUS"))
findings <- findings[order(findings$IID), ]
write.table(findings, FINDINGS, sep = "\t", quote = FALSE, row.names = FALSE)

is_high <- cmono %in% c("high_PTV", "high_missense")
etype <- ifelse(is_high, sample(c("DEE", "Focal"), n, TRUE, prob = c(0.85, 0.15)),
         ifelse(cmono == "incomplete", sample(c("Focal", "Generalised"), n, TRUE, prob = c(0.85, 0.15)),
                sample(c("Generalised", "Focal", "DEE"), n, TRUE, prob = c(0.45, 0.5, 0.05))))
onset <- round(ifelse(etype == "DEE", rgamma(n, 2, 2), rgamma(n, 4, 0.3)), 1)
sex <- sample(c("M", "F"), n, TRUE)
clin <- data.frame(`Participant ID` = iid,
                   Sex = ifelse(runif(n) < 0.1, ifelse(sex == "M", "male", "female"), sex),
                   `Age at onset (years)` = ifelse(runif(n) < 0.03, "unknown", as.character(onset)),
                   `Epilepsy type` = etype,
                   `Drug resistant` = ifelse(runif(n) < ifelse(etype == "DEE", 0.8, 0.3), "Yes", "No"),
                   check.names = FALSE)
clin$Sex[sample(n, 5)] <- tolower(clin$Sex[sample(n, 5)])
# 10 people with clinical data but no genotypes, and the duplicate sample has none
clin <- clin[clin$`Participant ID` != iid[n], ]
ghost <- data.frame(`Participant ID` = sprintf("EPI%04d", sample(setdiff(1000:9999, as.integer(sub("EPI", "", iid[-n]))), 10)),
                    Sex = sample(c("M", "F"), 10, TRUE), `Age at onset (years)` = "12",
                    `Epilepsy type` = "Focal", `Drug resistant` = "No", check.names = FALSE)
clin <- rbind(clin, ghost)
clin <- clin[sample(nrow(clin)), ]
write.table(clin, CLINICAL, sep = "\t", quote = FALSE, row.names = FALSE)

# 6. Base GWAS summary statistics
message("Writing GWAS summary statistics ...")
make_gwas <- function(eff_snp, eff_std, n_gwas) {
  n_gwas <- unname(n_gwas)
  marg <- numeric(M)
  for (k in seq_along(eff_snp)) {
    cols <- which(snp$block == snp$block[eff_snp[k]])
    r <- suppressWarnings(cor(geno_clean[, cols], geno_clean[, eff_snp[k]]))
    marg[cols] <- marg[cols] + ifelse(is.na(r), 0, r) * eff_std[k]
  }
  z <- marg * sqrt(n_gwas) + rnorm(M)
  f <- unname(p_pool)
  d <- data.frame(SNP = snp$id, CHR = snp$chr, BP = snp$bp, A1 = snp$A1, A2 = snp$A2,
                  FRQ = round(f, 4), BETA = signif(z / sqrt(n_gwas * 2 * f * (1 - f)), 4),
                  SE = signif(1 / sqrt(n_gwas * 2 * f * (1 - f)), 4), P = signif(2 * pnorm(-abs(z)), 4), N = n_gwas)
  flip <- runif(M) < 0.5                                  # effect allele is not always A1 in the cohort
  d[flip, c("A1", "A2")] <- d[flip, c("A2", "A1")]
  d$BETA[flip] <- -d$BETA[flip]; d$FRQ[flip] <- 1 - d$FRQ[flip]
  d <- d[-sample(M, round(0.03 * M)), ]                   # some cohort SNPs absent from the GWAS
  d <- rbind(d, d[sample(nrow(d), 10), ])                 # a few duplicated rows
  d$P[sample(nrow(d), 5)] <- NA                           # and a few missing P values
  d[order(d$CHR, d$BP), ]
}
gw_epi <- make_gwas(causal_snp, beta_std, N_GWAS["epilepsy"])
gw_ht  <- make_gwas(height_snp, rnorm(N_CAUSAL, 0, sqrt(0.5 / N_CAUSAL)), N_GWAS["height"])
gz <- function(d, f) { con <- gzfile(f, "w"); write.table(d, con, sep = "\t", quote = FALSE, row.names = FALSE); close(con) }
gz(gw_epi, GWAS[["epilepsy"]]); gz(gw_ht, GWAS[["height"]])

# 7. Truth for the workshop leader
truth <- data.frame(IID = iid, ancestry = c("A", "B")[cpop], rare_variant = cmono,
                    true_polygenic_liability = round(c(G[case], G[case][dup_from]), 3))
write.table(truth, file.path(DATA_DIR, "INSTRUCTOR_truth.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(qc_truth, file.path(DATA_DIR, "INSTRUCTOR_qc_problems.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

message("Done. Files in ", DATA_DIR, ":")
print(table(rare_variant = cmono, ancestry = c("A", "B")[cpop]))
message("True polygenic liability by rare-variant group (mean):")
print(round(tapply(truth$true_polygenic_liability, truth$rare_variant, mean), 2))
