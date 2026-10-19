#!/usr/bin/env Rscript
# Joint PRS model for the binary outcome with thresholds fixed from the case-control run: repeated CV and permutation tests.
# Usage: N_REP=100 N_PERM=1000 Rscript 07_cv_and_permutation.R   (about 5-10 min with the defaults)
suppressPackageStartupMessages({
    library(ggplot2)
})

WORKDIR   <- Sys.getenv("WORKDIR", "/re_gecip/<your_gecip>/<project>/<cohort>")
WORK      <- Sys.getenv("WITHIN_DIR", file.path(WORKDIR, "within_case"))
OUTCOME   <- Sys.getenv("OUTCOME", "outcome")
PT_SOURCE <- Sys.getenv("PT_SOURCE", file.path(WORKDIR, "results"))
RES_BIN   <- file.path(WORK, "results_binary")
EXTRA_COV <- Sys.getenv("EXTRA_COV_TSV", file.path(WORK, "extra_cov.txt"))
PFX       <- sprintf("within_case_%s", OUTCOME)

BASES    <- strsplit(Sys.getenv("JOINT_PRS_LIST", "epilepsy intelligence mvage bmi"), "[[:space:]]+")[[1]]
NB       <- length(BASES)
PTS      <- c(0.0001, 0.001, 0.01, 0.05, 0.1, 0.5, 1)     # bar levels used in every PRSice run
PT_FIXED <- NA          # e.g. 1 to use one threshold for every PRS instead of the case-control ones
N_REP    <- as.integer(Sys.getenv("N_REP",  "100"))
N_PERM   <- as.integer(Sys.getenv("N_PERM", "1000"))
N_FOLDS  <- 5
set.seed(42)

auc <- function(y, p) {                                    # Mann-Whitney AUC
    n1 <- sum(y == 1); n0 <- sum(y == 0)
    r  <- rank(c(p[y == 1], p[y == 0]))
    (sum(r[seq_len(n1)]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
emp_p <- function(obs, null) (1 + sum(null >= obs)) / (1 + length(null))
q95   <- function(x) quantile(x, c(0.5, 0.025, 0.975), names = FALSE)
say   <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), sprintf(...), "\n", sep = "")

cov <- read.table(file.path(WORK, "cov.txt"),          header = TRUE, stringsAsFactors = FALSE)
oc  <- read.table(file.path(WORK, "pheno_binary.txt"), header = TRUE, stringsAsFactors = FALSE)
oc$Y01 <- ifelse(oc$Outcome == 2, 1L, ifelse(oc$Outcome == 1, 0L, NA_integer_))
d <- merge(oc[!is.na(oc$Y01), c("IID", "Y01")], cov, by = "IID")

read_scores <- function(base) {                            # n x 7 matrix, one column per threshold
    x <- read.table(file.path(RES_BIN, sprintf("PRS_%s.all.score", base)),
                    header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
    if (!"IID" %in% names(x)) names(x)[1:2] <- c("FID", "IID")
    cols <- setdiff(names(x), c("FID", "IID"))
    vals <- as.numeric(sub("^Pt_", "", cols))
    m <- sapply(PTS, function(p) {
        i <- which.min(abs(vals - p))
        if (abs(vals[i] - p) > 1e-6 * p) stop(sprintf("%s: no Pt_%g column in .all.score", base, p))
        x[[cols[i]]]
    })
    rownames(m) <- x$IID
    m
}
raw <- lapply(setNames(BASES, BASES), read_scores)

n_all <- nrow(d)
d <- d[d$IID %in% Reduce(intersect, lapply(raw, rownames)) & d$Sex %in% c(1, 2), ]
if (nrow(d) < n_all) say("dropped %d samples without all PRSes or with unknown sex", n_all - nrow(d))
d$Male <- as.integer(d$Sex == 1)
y <- d$Y01
n <- nrow(d)
S <- lapply(raw, function(m) scale(m[d$IID, , drop = FALSE]))   # z-scored within this set
say("%s set: %d samples (%d yes, %d no)", OUTCOME, n, sum(y), n - sum(y))

PC_TERMS <- "PC1 + PC2 + PC3 + PC4"
F <- list(joint = as.formula(paste("Y01 ~", paste0("PRS_", BASES, collapse = " + "), "+", PC_TERMS, "+ Male")),
          base  = as.formula(paste("Y01 ~", PC_TERMS, "+ Male")),
          pcs   = as.formula(paste("Y01 ~", PC_TERMS)),
          sex   = Y01 ~ Male)

with_prs <- function(k, SS = S) {                          # PRS_<base> at threshold index k[base]
    dd <- d[, c("Y01", "PC1", "PC2", "PC3", "PC4", "Male")]
    for (b in BASES) dd[[paste0("PRS_", b)]] <- SS[[b]][, k[b]]
    dd
}
fitp <- function(f, dd, train = rep(TRUE, nrow(dd)), test = train) {
    fit <- suppressWarnings(glm(f, family = binomial, data = dd[train, ]))
    list(fit = fit, pred = suppressWarnings(predict(fit, newdata = dd[test, ], type = "response")))
}
single <- function(b, k, yy = y, SS = S, rows = seq_len(n)) {   # outcome ~ one PRS + PCs + sex
    dd <- data.frame(Y01 = yy[rows], PRS = SS[[b]][rows, k], d[rows, c("PC1", "PC2", "PC3", "PC4", "Male")])
    coef(summary(suppressWarnings(glm(Y01 ~ ., family = binomial, data = dd))))["PRS", ]
}
# PRSice's rule: the threshold where the PRS is most significant for the outcome, given PCs and sex
pick_pt <- function(rows = seq_len(n), yy = y, SS = S)
    sapply(BASES, function(b) which.min(sapply(seq_along(PTS), function(k) single(b, k, yy, SS, rows)[[4]])))

# Thresholds fixed in advance (case-control run) are the main analysis. Picking them
# on the outcome in these same patients is kept as a sensitivity analysis, with the
# pick redone inside every CV fold and every permutation.
pt_fixed <- if (is.na(PT_FIXED)) {
    sapply(BASES, function(b) read.table(file.path(PT_SOURCE, sprintf("PRS_%s.summary", b)),
                                         header = TRUE)$Threshold[1])
} else setNames(rep(PT_FIXED, length(BASES)), BASES)
k_fix  <- sapply(pt_fixed, function(p) which.min(abs(PTS - p)))
k_pick <- pick_pt()
pt_prsice <- sapply(BASES, function(b) read.table(file.path(RES_BIN, sprintf("PRS_%s.summary", b)),
                                                  header = TRUE)$Threshold[1])
nsnp_at <- function(b, k) {
    f <- file.path(RES_BIN, sprintf("PRS_%s.prsice", b))
    if (!file.exists(f)) return(NA)
    x <- read.table(f, header = TRUE)
    x$Num_SNP[which.min(abs(x$Threshold - PTS[k]))]
}
say("thresholds fixed in advance: %s", paste(sprintf("%s %g", BASES, PTS[k_fix]), collapse = ", "))
say("thresholds picked on %s   : %s", OUTCOME, paste(sprintf("%s %g", BASES, PTS[k_pick]), collapse = ", "))
if (any(abs(PTS[k_pick] - pt_prsice) > 1e-9))
    say("note: PRSice's own pick differs for %s (it may have used other covariates)",
        paste(BASES[abs(PTS[k_pick] - pt_prsice) > 1e-9], collapse = ", "))

D_fix  <- with_prs(k_fix)
D_pick <- with_prs(k_pick)
fb     <- fitp(F$base, D_fix)                               # the baseline does not involve the PRSes
insample <- c(pcs  = auc(y, fitp(F$pcs, D_fix)$pred),
              base = auc(y, fb$pred),
              sex  = auc(y, fitp(F$sex, D_fix)$pred),
              joint_fixed  = auc(y, fitp(F$joint, D_fix)$pred),
              joint_picked = auc(y, fitp(F$joint, D_pick)$pred))
say("in-sample AUC, thresholds picked on %s: %.3f", OUTCOME, insample["joint_picked"])

lr_gain <- function(k, SS = S) {
    fj <- fitp(F$joint, with_prs(k, SS))
    c(LR = as.numeric(2 * (logLik(fj$fit) - logLik(fb$fit))), dAUC = auc(y, fj$pred) - auc(y, fb$pred))
}
obs_fixed  <- lr_gain(k_fix)
obs_picked <- lr_gain(k_pick)

# Repeated stratified 5-fold CV; all models share the same folds in each repeat
make_folds <- function(yy) {
    f <- integer(length(yy))
    for (cl in 0:1) { i <- which(yy == cl); f[i] <- sample(rep_len(seq_len(N_FOLDS), length(i))) }
    f
}
cv_rep <- function(folds) {
    pred <- matrix(NA_real_, n, 5, dimnames = list(NULL, names(insample)))
    for (f in seq_len(N_FOLDS)) {
        tr <- folds != f; te <- folds == f
        pred[te, "pcs"]  <- fitp(F$pcs,  D_fix, tr, te)$pred
        pred[te, "base"] <- fitp(F$base, D_fix, tr, te)$pred
        pred[te, "sex"]  <- fitp(F$sex,  D_fix, tr, te)$pred
        pred[te, "joint_fixed"]  <- fitp(F$joint, D_fix, tr, te)$pred
        pred[te, "joint_picked"] <- fitp(F$joint, with_prs(pick_pt(which(tr))), tr, te)$pred   # re-picked on the training folds
    }
    apply(pred, 2, function(p) auc(y, p))
}
say("cross-validation: %d repeats x %d folds", N_REP, N_FOLDS)
cv <- matrix(NA_real_, N_REP, 5, dimnames = list(NULL, names(insample)))
for (r in seq_len(N_REP)) {
    cv[r, ] <- cv_rep(make_folds(y))
    if (r %% 20 == 0) say("  %d / %d", r, N_REP)
}

# Outcome labels shuffled: null for the baselines and for the per-PRS threshold search
say("permutations of %s labels: %d", OUTCOME, N_PERM)
null_auc  <- matrix(NA_real_, N_PERM, 3, dimnames = list(NULL, c("pcs", "base", "sex")))
null_minp <- matrix(NA_real_, N_PERM, length(BASES), dimnames = list(NULL, BASES))
for (i in seq_len(N_PERM)) {
    yp <- sample(y)
    Dp <- D_fix; Dp$Y01 <- yp
    for (m in colnames(null_auc)) null_auc[i, m] <- auc(yp, fitp(F[[m]], Dp)$pred)
    for (b in BASES) null_minp[i, b] <- min(sapply(seq_along(PTS), function(k) single(b, k, yp)[[4]]))
    if (i %% 200 == 0) say("  %d / %d", i, N_PERM)
}

# The PRS residuals (after PCs + sex) are shuffled together across patients, so each
# shuffled PRS keeps its link with ancestry and sex but loses any link with the outcome.
X  <- qr(cbind(1, as.matrix(d[, c("PC1", "PC2", "PC3", "PC4", "Male")])))
FP <- lapply(S, function(m) qr.fitted(X, m))
RP <- lapply(S, function(m) qr.resid(X, m))
say("permutations of the PRS block: %d", N_PERM)
null_fixed <- null_picked <- matrix(NA_real_, N_PERM, 2, dimnames = list(NULL, c("LR", "dAUC")))
for (i in seq_len(N_PERM)) {
    o   <- sample(n)
    SSp <- mapply(function(fp, rp) fp + rp[o, , drop = FALSE], FP, RP, SIMPLIFY = FALSE)
    null_fixed[i, ]  <- lr_gain(k_fix, SSp)
    null_picked[i, ] <- lr_gain(pick_pt(yy = y, SS = SSp), SSp)   # threshold search redone every time
    if (i %% 200 == 0) say("  %d / %d", i, N_PERM)
}

out <- function(x, name) {
    write.table(x, file.path(WORK, name), sep = "\t", row.names = FALSE, quote = FALSE)
    say("wrote %s", name)
}
sig <- function(x) signif(x, 3)

thr <- do.call(rbind, lapply(BASES, function(b) {
    fx <- single(b, k_fix[b])
    obs_minp <- min(sapply(seq_along(PTS), function(k) single(b, k)[[4]]))
    data.frame(PRS = b,
               Pt_fixed = PTS[k_fix[b]], N_SNP_fixed = nsnp_at(b, k_fix[b]),
               OR_perSD_fixed = sig(exp(fx[[1]])), P_fixed = sig(fx[[4]]),
               Pt_picked_on_outcome = PTS[k_pick[b]], Pt_PRSice = pt_prsice[b],
               P_picked_nominal = sig(obs_minp),
               P_picked_corrected = sig((1 + sum(null_minp[, b] <= obs_minp)) / (1 + N_PERM)))
}))
out(thr, paste0(PFX, "_thresholds.tsv"))

co <- coef(summary(fitp(F$joint, D_fix)$fit))
jc <- do.call(rbind, lapply(BASES, function(b) {
    r <- co[paste0("PRS_", b), ]
    data.frame(PRS = b, Pt = PTS[k_fix[b]], N_SNP = nsnp_at(b, k_fix[b]),
               beta_z = sig(r[1]), SE = sig(r[2]), OR_perSD = sig(exp(r[1])),
               OR_2.5 = sig(exp(r[1] - 1.96 * r[2])), OR_97.5 = sig(exp(r[1] + 1.96 * r[2])), P = sig(r[4]))
}))
out(jc, paste0(PFX, "_joint_coef.tsv"))

LABELS <- c(pcs = "PCs only", base = "PCs + sex", sex = "Sex only",
            joint_fixed  = sprintf("%d PRSes + PCs + sex (thresholds fixed in advance)", NB),
            joint_picked = sprintf("%d PRSes + PCs + sex (thresholds picked on %s)", NB, OUTCOME))
perm_p <- c(pcs  = emp_p(insample["pcs"],  null_auc[, "pcs"]),
            base = emp_p(insample["base"], null_auc[, "base"]),
            sex  = emp_p(insample["sex"],  null_auc[, "sex"]),
            joint_fixed  = emp_p(obs_fixed["LR"],  null_fixed[, "LR"]),
            joint_picked = emp_p(obs_picked["LR"], null_picked[, "LR"]))
models <- do.call(rbind, lapply(names(LABELS), function(m) {
    s <- q95(cv[, m])
    data.frame(Model = LABELS[[m]], AUC_insample = sig(insample[m]),
               AUC_CV_median = sig(s[1]), AUC_CV_2.5 = sig(s[2]), AUC_CV_97.5 = sig(s[3]),
               Optimism = sig(insample[m] - s[1]),
               Permutation_test = if (startsWith(m, "joint")) "PRSes shuffled (added value over PCs + sex)" else "outcome labels shuffled",
               Perm_P = sig(perm_p[m]))
}))
out(models, paste0(PFX, "_models.tsv"))

inc <- do.call(rbind, lapply(c("fixed", "picked"), function(v) {
    obs  <- if (v == "fixed") obs_fixed else obs_picked
    nul  <- if (v == "fixed") null_fixed else null_picked
    gain <- cv[, paste0("joint_", v)] - cv[, "base"]
    s    <- q95(gain)
    data.frame(Thresholds = if (v == "fixed") "fixed in advance (case-control run)" else sprintf("picked on %s (search redone in every permutation and fold)", OUTCOME),
               LR_chisq = sig(obs["LR"]), df = NB,
               LR_P_parametric = if (v == "fixed") sig(pchisq(obs["LR"], NB, lower.tail = FALSE)) else NA,
               LR_P_permutation = sig(emp_p(obs["LR"], nul[, "LR"])),
               AUC_gain_insample = sig(obs["dAUC"]),
               AUC_gain_CV_median = sig(s[1]), AUC_gain_CV_2.5 = sig(s[2]), AUC_gain_CV_97.5 = sig(s[3]),
               Pct_CV_repeats_PRS_help = round(100 * mean(gain > 0)))
}))
out(inc, paste0(PFX, "_prs_increment.tsv"))

# Which PCs carry the baseline signal
Dz <- D_fix
for (pc in paste0("PC", 1:4)) Dz[[pc]] <- as.numeric(scale(Dz[[pc]]))
cp <- coef(summary(glm(F$pcs, family = binomial, data = Dz)))
pcs_tab <- do.call(rbind, lapply(paste0("PC", 1:4), function(pc) {
    a <- auc(y, Dz[[pc]])
    data.frame(PC = pc, OR_perSD = sig(exp(cp[pc, 1])),
               OR_2.5 = sig(exp(cp[pc, 1] - 1.96 * cp[pc, 2])), OR_97.5 = sig(exp(cp[pc, 1] + 1.96 * cp[pc, 2])),
               P = sig(cp[pc, 4]), AUC_alone = sig(max(a, 1 - a)))
}))
out(pcs_tab, paste0(PFX, "_pcs.tsv"))

# Numeric extra covariates are used as measurements; give categories as text labels.
if (file.exists(EXTRA_COV)) {
    ex <- read.table(EXTRA_COV, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
    ex <- ex[match(d$IID, ex$IID), setdiff(names(ex), "IID"), drop = FALSE]
    ex <- ex[, sapply(ex, function(v) length(unique(na.omit(v))) > 1), drop = FALSE]   # drop constant columns
    ok <- complete.cases(ex)
    say("extra covariates: %s (complete for %d of %d)", paste(names(ex), collapse = ", "), sum(ok), n)

    assoc <- do.call(rbind, lapply(names(ex), function(v) do.call(rbind, lapply(paste0("PC", 1:4), function(pc) {
        x <- ex[ok, v]; p <- d[ok, pc]
        if (is.numeric(x)) {
            t <- cor.test(p, x)
            data.frame(Covariate = v, PC = pc, Measure = "Pearson r", Value = sig(t$estimate), P = sig(t$p.value))
        } else {
            a <- anova(lm(p ~ factor(x)))
            data.frame(Covariate = v, PC = pc, Measure = "R2", Value = sig(a[1, 2] / sum(a[, 2])), P = sig(a[1, 5]))
        }
    }))))
    out(assoc, paste0(PFX, "_pc_vs_extra.tsv"))

    XE <- model.matrix(~ ., data = ex[ok, , drop = FALSE])[, -1, drop = FALSE]
    colnames(XE) <- make.names(colnames(XE))
    De <- cbind(D_fix[ok, c("Y01", "PC1", "PC2", "PC3", "PC4", "Male")], XE)
    ye <- y[ok]
    FE <- list(extra     = as.formula(paste("Y01 ~ Male +", paste(colnames(XE), collapse = " + "))),
               extra_pcs = as.formula(paste("Y01 ~ Male +", paste(colnames(XE), collapse = " + "), "+", PC_TERMS)),
               pcs_sex   = F$base)
    fe  <- fitp(FE$extra, De); fep <- fitp(FE$extra_pcs, De)
    lr  <- as.numeric(2 * (logLik(fep$fit) - logLik(fe$fit)))
    ins <- sapply(FE, function(f) auc(ye, fitp(f, De)$pred))
    cve <- t(replicate(N_REP, {
        folds <- make_folds(ye)
        pred  <- matrix(NA_real_, sum(ok), 3, dimnames = list(NULL, names(FE)))
        for (f in seq_len(N_FOLDS)) for (m in names(FE))
            pred[folds == f, m] <- fitp(FE[[m]], De, folds != f, folds == f)$pred
        apply(pred, 2, function(p) auc(ye, p))
    }))
    ext <- do.call(rbind, lapply(names(FE), function(m) {
        s <- q95(cve[, m])
        data.frame(Model = c(extra = "extra covariates + sex", extra_pcs = "extra covariates + sex + PCs",
                             pcs_sex = "PCs + sex")[[m]],
                   N = sum(ok), AUC_insample = sig(ins[m]),
                   AUC_CV_median = sig(s[1]), AUC_CV_2.5 = sig(s[2]), AUC_CV_97.5 = sig(s[3]),
                   PCs_LR_given_extra = if (m == "extra_pcs") sig(lr) else NA,
                   PCs_P_given_extra  = if (m == "extra_pcs") sig(pchisq(lr, 4, lower.tail = FALSE)) else NA)
    }))
    out(ext, paste0(PFX, "_extra_cov.tsv"))

    # PRS block over a baseline that includes the extra covariates: thresholds fixed in
    # advance, same residual permutation and repeated CV as above.
    PRS_TERMS <- paste0("PRS_", BASES, collapse = " + ")
    for (b in BASES) De[[paste0("PRS_", b)]] <- D_fix[ok, paste0("PRS_", b)]
    ne <- sum(ok)
    BL <- c(extra     = paste("Y01 ~ Male +", paste(colnames(XE), collapse = " + ")),
            extra_pcs = paste("Y01 ~ Male +", paste(colnames(XE), collapse = " + "), "+", PC_TERMS))
    BL_LAB <- c(extra = "extra covariates + sex", extra_pcs = "extra covariates + sex + PCs")
    lr_e <- function(fbm, fj_f, DD) {
        fj <- fitp(fj_f, DD)
        c(LR = as.numeric(2 * (logLik(fj$fit) - logLik(fbm$fit))), dAUC = auc(ye, fj$pred) - auc(ye, fbm$pred))
    }
    say("PRS block given the extra covariates: %d permutations and %d CV repeats per baseline", N_PERM, N_REP)
    pge <- list(); pgc <- list()
    for (m in names(BL)) {
        fb_f <- as.formula(BL[[m]])
        fj_f <- as.formula(paste(BL[[m]], "+", PRS_TERMS))
        fbm  <- fitp(fb_f, De)
        obs  <- lr_e(fbm, fj_f, De)
        Xb   <- qr(model.matrix(fb_f, De))
        Pm   <- sapply(BASES, function(b) De[[paste0("PRS_", b)]])
        FPm  <- qr.fitted(Xb, Pm); RPm <- qr.resid(Xb, Pm)
        nul  <- t(replicate(N_PERM, {
            Dp <- De
            Dp[paste0("PRS_", BASES)] <- as.data.frame(FPm + RPm[sample(ne), , drop = FALSE])
            lr_e(fbm, fj_f, Dp)
        }))
        cvg <- t(replicate(N_REP, {
            folds <- make_folds(ye)
            pb <- pj <- numeric(ne)
            for (f in seq_len(N_FOLDS)) {
                tr <- folds != f; te <- folds == f
                pb[te] <- fitp(fb_f, De, tr, te)$pred
                pj[te] <- fitp(fj_f, De, tr, te)$pred
            }
            c(base = auc(ye, pb), joint = auc(ye, pj))
        }))
        gain_e <- cvg[, "joint"] - cvg[, "base"]
        s <- q95(gain_e)
        pge[[m]] <- data.frame(Baseline = BL_LAB[[m]], N = ne, Cases = sum(ye),
                               LR_chisq = sig(obs["LR"]), df = length(BASES),
                               LR_P_parametric = sig(pchisq(obs["LR"], length(BASES), lower.tail = FALSE)),
                               LR_P_permutation = sig(emp_p(obs["LR"], nul[, "LR"])),
                               AUC_CV_baseline = sig(median(cvg[, "base"])), AUC_CV_with_PRS = sig(median(cvg[, "joint"])),
                               AUC_gain_insample = sig(obs["dAUC"]),
                               AUC_gain_CV_median = sig(s[1]), AUC_gain_CV_2.5 = sig(s[2]), AUC_gain_CV_97.5 = sig(s[3]),
                               Pct_CV_repeats_PRS_help = round(100 * mean(gain_e > 0)))
        ce <- coef(summary(fitp(fj_f, De)$fit))
        pgc[[m]] <- do.call(rbind, lapply(BASES, function(b) {
            r <- ce[paste0("PRS_", b), ]
            data.frame(Baseline = BL_LAB[[m]], PRS = b, Pt = PTS[k_fix[b]], OR_perSD = sig(exp(r[1])),
                       OR_2.5 = sig(exp(r[1] - 1.96 * r[2])), OR_97.5 = sig(exp(r[1] + 1.96 * r[2])), P = sig(r[4]))
        }))
    }
    pge <- do.call(rbind, pge); out(pge, paste0(PFX, "_prs_given_extra.tsv"))
    pgc <- do.call(rbind, pgc); out(pgc, paste0(PFX, "_prs_given_extra_coef.tsv"))
} else {
    say("no %s: skipping the extra-covariate analyses", EXTRA_COV)
}

ORDER <- c("joint_picked", "joint_fixed", "sex", "base", "pcs")      # bottom to top
YLAB  <- c(pcs = "PCs only", base = "PCs + sex", sex = "Sex only",
           joint_fixed  = sprintf("%d PRSes + PCs + sex\n(thresholds fixed in advance)", NB),
           joint_picked = sprintf("%d PRSes + PCs + sex\n(thresholds picked on %s)", NB, OUTCOME))
cv_long <- data.frame(model = factor(rep(colnames(cv), each = N_REP), levels = ORDER), AUC = as.vector(cv),
                      what = sprintf("Cross-validated AUC (%d repeats of 5-fold)", N_REP))
ins_df  <- data.frame(model = factor(names(insample), levels = ORDER), AUC = insample,
                      what = "In-sample AUC (fit and scored on the same patients)")

g <- ggplot(cv_long, aes(AUC, model)) +
    geom_vline(xintercept = 0.5, linetype = "dashed", colour = "grey65", linewidth = 0.4) +
    geom_boxplot(aes(fill = what), colour = "#2F6FDB", width = 0.5, linewidth = 0.5,
                 outlier.shape = NA, alpha = 0.25) +
    geom_point(data = ins_df, aes(shape = what), colour = "#D1495B", fill = "white", size = 3.2, stroke = 1.1) +
    scale_fill_manual(values = "#2F6FDB", name = NULL) +
    scale_shape_manual(values = 23, name = NULL) +
    scale_y_discrete(labels = YLAB) +
    labs(title = sprintf("Within-case %s: how well each model predicts patients it has not seen (%d yes, %d no)",
                         OUTCOME, sum(y), n - sum(y)),
         subtitle = "Dashed line = chance (AUC 0.5). The gap between a diamond and its box is the optimism of scoring in-sample.",
         x = "AUC", y = NULL) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 12.5), plot.title.position = "plot",
          plot.subtitle = element_text(size = 9.5, colour = "grey30", margin = margin(b = 8)),
          legend.position = "bottom", legend.box = "horizontal",
          legend.text = element_text(size = 9.5), panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank(), axis.text.y = element_text(size = 10.5, colour = "grey15"))

ggsave(file.path(WORK, paste0(PFX, "_cv_auc.png")), g, width = 9.5, height = 5.8, dpi = 300)
ggsave(file.path(WORK, paste0(PFX, "_cv_auc.pdf")), g, width = 9.5, height = 5.8)
say("wrote %s_cv_auc.{png,pdf}", PFX)

cat("\nModels:\n");                         print(models, row.names = FALSE)
cat("\nDo the PRSes add anything?\n");      print(inc, row.names = FALSE)
cat("\nThresholds:\n");                     print(thr, row.names = FALSE)
cat("\nPCs:\n");                            print(pcs_tab, row.names = FALSE)
if (exists("pge")) {
    cat("\nPRSes over the extra covariates:\n"); print(pge, row.names = FALSE)
    cat("\nPer-PRS ORs, adjusted models:\n");    print(pgc, row.names = FALSE)
}
