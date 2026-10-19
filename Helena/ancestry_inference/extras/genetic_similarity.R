#!/usr/bin/env Rscript
# Proportional genetic similarity to 21 reference populations (Prive et al. 2022) from a PLINK bed set.
# Usage: Rscript genetic_similarity.R <bed_prefix> <ref_dir with ref_freqs.csv.gz and projection.csv.gz> [out_dir]

library(bigsnpr)
library(bigreadr)
library(quadprog)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript genetic_similarity.R <bed_prefix> <ref_dir> [out_dir]")
}
bed_prefix <- args[1]
REF_DIR    <- args[2]
OUT_DIR    <- if (length(args) >= 3) args[3] else "gen_sim"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("[", format(Sys.time()), "] Input:", bed_prefix, "\n")

all_freq <- fread2(file.path(REF_DIR, "ref_freqs.csv.gz"))
cat("Reference: ", nrow(all_freq), " variants x ",
    ncol(all_freq) - 5, " populations\n", sep = "")

projection <- fread2(file.path(REF_DIR, "projection.csv.gz"))
cat("Loadings: ", nrow(projection), " variants x ",
    ncol(projection) - 5, " PCs\n", sep = "")

pop_names <- colnames(all_freq)[-(1:5)]

bim <- fread2(paste0(bed_prefix, ".bim"),
              select = c(1, 4, 5, 6),
              col.names = c("chr", "pos", "a1", "a0"))
bim$chr <- as.integer(gsub("chr", "", bim$chr))
bim$beta <- 1
cat("Study variants: ", nrow(bim), "\n", sep = "")

# beta comes back as -1 for variants matched with alleles swapped, which flips
# the centring and scaling in the projection below.
matched <- snp_match(bim, all_freq[, 1:5])
cat("Matched variants: ", nrow(matched), "\n", sep = "")

bed_file <- paste0(bed_prefix, ".bed")
rds_file <- snp_readBed2(bed_file, ind.col = matched$`_NUM_ID_.ss`,
                         backingfile = file.path(OUT_DIR, "study_matched"))
obj <- snp_attach(rds_file)
G <- obj$genotypes
cat("Genotypes: ", nrow(G), " samples x ", ncol(G), " variants\n", sep = "")

nb_na <- big_counts(G)[4, ]
miss_rate <- nb_na / nrow(G)
ind_keep <- which(miss_rate < 0.05)
cat("Variants with < 5% missing: ", length(ind_keep), " / ", ncol(G), "\n", sep = "")

G2 <- snp_fastImputeSimple(G)

# Shrinkage correction for projected PCs 1-16, from Prive et al. 2022
correction <- c(1, 1, 1, 1.008, 1.021, 1.034, 1.052, 1.074,
                1.099, 1.123, 1.15, 1.195, 1.256, 1.321, 1.382, 1.443)

PROJ <- as.matrix(projection[matched$`_NUM_ID_`[ind_keep], -(1:5)])

all_proj <- big_prodMat(
  G2,
  sweep(PROJ, 2, correction / 2, '*'),
  ind.col = ind_keep,
  center = 1 - matched$beta[ind_keep],
  scale = matched$beta[ind_keep]
)

# Population centres in the same PC space (16 PCs x 21 populations)
X <- crossprod(
  as.matrix(projection[matched$`_NUM_ID_`[ind_keep], -(1:5)]),
  as.matrix(all_freq[matched$`_NUM_ID_`[ind_keep], -(1:5)])
)

group <- pop_names
grp_fct <- factor(group, unique(group))

# Non-negative proportions summing to 1, by quadratic programming
cp_X_pd <- Matrix::nearPD(crossprod(X), base.matrix = TRUE)
Amat <- cbind(1, diag(ncol(X)))
bvec <- c(1, rep(0, ncol(X)))

all_res <- apply(all_proj, 1, function(y) {
  sol <- quadprog::solve.QP(
    Dmat = cp_X_pd$mat,
    dvec = crossprod(y, X),
    Amat = Amat,
    bvec = bvec,
    meq = 1
  )$sol
  tapply(sol, grp_fct, sum) |> round(7)
})

sim_props <- as.data.frame(t(all_res))

all_centres <- t(X)

all_sq_dist <- apply(all_centres, 1, function(centre) {
  rowSums(sweep(all_proj, 2, centre, '-')^2)
})
colnames(all_sq_dist) <- pop_names

# Squared PC distance rescaled to approximate FST (largest between-centre distance = 0.16)
max_dist <- max(dist(all_centres)^2)
fst_approx <- all_sq_dist * 0.16 / max_dist

thr_fst <- 0.002
thr_sq_dist <- max_dist * thr_fst / 0.16

closest_group <- apply(all_sq_dist, 1, function(x) {
  idx <- which.min(x)
  if (x[idx] < thr_sq_dist) pop_names[idx] else NA_character_
})

fam <- fread2(paste0(bed_prefix, ".fam"),
              col.names = c("FID", "IID", "PAT", "MAT", "SEX", "PHENO"))

out_sim <- data.frame(
  FID = fam$FID,
  IID = fam$IID,
  sim_props,
  closest_group = closest_group,
  stringsAsFactors = FALSE
)

fst_df <- as.data.frame(fst_approx)
colnames(fst_df) <- paste0("fst_", colnames(fst_df))
out_full <- cbind(out_sim, fst_df)

write.csv(out_full,
          file.path(OUT_DIR, "genetic_similarity_results.csv"),
          row.names = FALSE)

cat("Samples: ", nrow(out_full), "\n")
cat("Variants used: ", length(ind_keep), "\n")
cat("Closest group counts:\n")
print(table(closest_group, exclude = NULL))

# make.names() turns the parentheses in population names into double dots
eur_cols <- c("United.Kingdom", "Ireland", "Scandinavia", "Finland",
              "Italy", "Europe..South.East.", "Europe..South.West.",
              "Europe..North.East.", "Ashkenazi")
eur_cols_found <- intersect(names(sim_props), eur_cols)
if (length(eur_cols_found) > 0) {
  total_eur <- rowSums(sim_props[, eur_cols_found, drop = FALSE])
  cat("Mean European proportion: ", round(mean(total_eur) * 100, 1), "%\n", sep = "")
  cat("Samples with >80% / >70% / >50% European: ",
      sum(total_eur > 0.80), " / ", sum(total_eur > 0.70), " / ", sum(total_eur > 0.50),
      "\n", sep = "")

  eur_summary <- data.frame(
    FID = fam$FID,
    IID = fam$IID,
    total_eur_proportion = round(total_eur, 4),
    closest_group = closest_group,
    stringsAsFactors = FALSE
  )
  eur_summary <- eur_summary[order(-eur_summary$total_eur_proportion), ]
  write.csv(eur_summary,
            file.path(OUT_DIR, "eur_proportion_summary.csv"),
            row.names = FALSE)
}

cat("[", format(Sys.time()), "] Done, output in ", OUT_DIR, "\n", sep = "")
