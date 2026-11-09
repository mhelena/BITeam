# Helper functions used by the notebooks.
# You do not need to edit this file, but feel free to read it.
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2); library(stringr); library(patchwork)
})

# Colours used in every plot, so a group keeps its colour throughout
GROUP_COLOURS <- c("No variant reported" = "#7F7F7F", "VUS only" = "#E6A700", "Monogenic" = "#B2182B",
                   "High penetrance" = "#B2182B", "Incomplete penetrance" = "#2E5FA8", "Other gene" = "#7F7F7F",
                   "PTV" = "#B2182B", "Missense" = "#3E8E7E", "Case group" = "#B2182B", "Comparison group" = "#2E5FA8")

theme_workshop <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"),
          plot.subtitle = element_text(colour = "grey35"), legend.position = "bottom")
}

# Running plink2 and PRSice from R
# Runs a command-line tool, keeps everything it prints in results/logs/<log_name>.txt,
# and stops with a readable message if the tool fails.
run_tool <- function(tool, args, log_name, show = 8) {
  if (Sys.which(tool) == "" && !file.exists(tool))
    stop("Cannot find '", tool, "'. Open config.R and set the full path to it.", call. = FALSE)
  dir.create(file.path(RESULTS, "logs"), showWarnings = FALSE, recursive = TRUE)
  log <- file.path(RESULTS, "logs", paste0(log_name, ".txt"))
  status <- suppressWarnings(system2(tool, args, stdout = log, stderr = log))
  out <- readLines(log, warn = FALSE)
  out <- out[!grepl("Progress|Processing [0-9]|^\\s*$|\r", out)]           # drop progress counters
  cat(paste0("$ ", basename(tool), " ", paste(args, collapse = " ")), "\n...\n")
  cat(tail(out, show), sep = "\n")
  if (!identical(as.integer(status), 0L))
    stop(basename(tool), " failed (exit code ", status, "). Read the full log: ", log, call. = FALSE)
  invisible(log)
}

# plink2 writes headers that start with '#'; read them as normal column names
read_plink <- function(file) {
  d <- read_tsv(file, show_col_types = FALSE, comment = "")
  names(d) <- sub("^#", "", names(d))
  d
}

# Comparing genotype groups inside the cohort
# case / control: conditions written with column names, as in dplyr::filter().
# Example:
#   compare_groups(dat, case = dx_status == "No variant reported", control = dx_status == "Monogenic")
# Model: PRS (in SD) ~ group + covariates. 'diff_SD' is how many standard deviations
# higher the case group's PRS is than the comparison group's, after adjustment.
compare_groups <- function(data, case, control, prs = "epilepsy_z",
                           covariates = c("sex", paste0("PC", 1:N_PCS)), name = NULL) {
  case <- rlang::enquo(case); control <- rlang::enquo(control)
  d <- data %>%
    mutate(.case = coalesce(!!case, FALSE), .control = coalesce(!!control, FALSE))
  both <- sum(d$.case & d$.control)
  if (both > 0) stop(both, " people meet both the case and the comparison definitions.", call. = FALSE)
  d <- d %>% filter(.case | .control) %>% mutate(group = as.integer(.case))
  if (min(table(factor(d$group, 0:1))) < 3) stop("Fewer than 3 people in one of the groups.", call. = FALSE)
  fit <- lm(reformulate(c("group", covariates), prs), data = d)
  est <- summary(fit)$coefficients["group", ]
  ci <- confint(fit)["group", ]
  tibble(contrast  = if (is.null(name)) paste(rlang::as_label(case), "vs", rlang::as_label(control)) else name,
         prs       = prs,
         n_case    = sum(d$group == 1),
         n_control = sum(d$group == 0),
         diff_SD   = unname(est["Estimate"]),
         low       = unname(ci[1]),
         high      = unname(ci[2]),
         P         = unname(est["Pr(>|t|)"]))
}

# Violin + box + points for one comparison, PRS adjusted for the covariates
plot_groups <- function(data, case, control, prs = "epilepsy_z",
                        covariates = c("sex", paste0("PC", 1:N_PCS)), title = NULL, labels = NULL) {
  case <- rlang::enquo(case); control <- rlang::enquo(control)
  d <- data %>%
    mutate(.case = coalesce(!!case, FALSE), .control = coalesce(!!control, FALSE)) %>%
    filter(.case | .control)
  d$adj <- resid(lm(reformulate(covariates, prs), data = d)) + mean(d[[prs]])
  d$grp <- factor(ifelse(d$.case, "Case group", "Comparison group"), c("Comparison group", "Case group"))
  if (is.null(labels)) labels <- c(rlang::as_label(control), rlang::as_label(case))
  lab <- c(paste0(labels[1], "\n(n = ", sum(!d$.case), ")"), paste0(labels[2], "\n(n = ", sum(d$.case), ")"))
  ggplot(d, aes(grp, adj, fill = grp, colour = grp)) +
    geom_violin(alpha = 0.15, colour = NA) +
    geom_jitter(width = 0.12, size = 1.1, alpha = 0.45) +
    geom_boxplot(width = 0.18, alpha = 0.6, colour = "black", outlier.shape = NA) +
    scale_x_discrete(labels = lab) +
    scale_fill_manual(values = GROUP_COLOURS, guide = "none") +
    scale_colour_manual(values = GROUP_COLOURS, guide = "none") +
    labs(x = NULL, y = paste(prs, "(SD, adjusted)"), title = title) +
    theme_workshop()
}

# Forest plot of several comparisons (the output of compare_groups, row-bound)
forest_plot <- function(res, title = "PRS differences between genotype groups") {
  res %>%
    mutate(label = paste0(contrast, "  (", n_case, " vs ", n_control, ")"),
           label = factor(label, rev(unique(label)))) %>%
    ggplot(aes(diff_SD, label, colour = prs)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_errorbarh(aes(xmin = low, xmax = high), height = 0.2, position = position_dodge(width = 0.5)) +
    geom_point(size = 2.8, position = position_dodge(width = 0.5)) +
    labs(x = "Difference in PRS, case group minus comparison group (SD, 95% CI)", y = NULL,
         colour = NULL, title = title) +
    theme_workshop()
}

# Name of the standardised score column for a trait at a P value threshold, e.g. score_col(dat, "epilepsy", 0.1)
score_col <- function(data, trait, pt) {
  cols <- grep(paste0("^", trait, "_Pt_.*_z$"), names(data), value = TRUE)
  pts <- as.numeric(sub(paste0("^", trait, "_Pt_(.*)_z$"), "\\1", cols))
  hit <- cols[abs(pts / pt - 1) < 1e-6]
  if (length(hit) != 1) stop("No ", trait, " score at threshold ", pt, call. = FALSE)
  hit
}
