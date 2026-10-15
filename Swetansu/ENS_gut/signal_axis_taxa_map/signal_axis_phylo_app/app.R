# =============================================================================
# signal_axis_phylo_app.R
#
# Interactive highlighting of signal axes on the microbial phylogeny drawn by
# microbe_phylogeny.py (signal_axis_map_phylo).
#
# Run:   shiny::runApp("signal_axis_phylo_app.R")
#        (to browse the 16S-fallback tree instead: Sys.setenv(SIGNAL_AXIS_PREFIX = "signal_axis_map_phylo_16S_fallback_test") first)
# Needs: shiny, ggplot2, dplyr, readr, tidyr   (no ape / ggtree)
#
# Data (put these three files, plus the optional fourth, in a `data/` folder next
# to this script, or in the same folder):
#   signal_axis_map_phylo_tips.csv        one row per tip: label, GTDB species, rank, flags
#   taxa_axis_membership.csv              taxon x axis table: axis index/name/direction/colour
#   signal_axis_map_phylo.nwk             the pruned tree (tip labels = source_id)
#   signal_axis_map_phylo_provenance.json optional: route/release, marks arbitrary stems
#
# The tips CSV carries no axis membership and no tree, so the app joins it to the
# membership table on `taxid` and to the tree on `source_id`.
#
# What "breadth" means: the total length of the smallest subtree connecting the
# selected taxa (unrooted), in the tree's own branch-length units. It is NOT Faith's
# PD (which adds the path to the root). GTDB's bacterial and archaeal trees are
# inferred separately and joined here by arbitrary stems, so breadth is computed
# WITHIN each of those trees and summed: a lone taxon in one tree (the archaeon)
# contributes 0, and the arbitrary stems never count. A single taxon has breadth 0.
#
# The taxa are representative examples named in panel_curation.py, not taxa
# measured in any cohort: nothing here is an abundance or prevalence result.
# =============================================================================
suppressPackageStartupMessages({
  library(shiny); library(ggplot2); library(dplyr); library(readr); library(tidyr)
})

PREFIX <- Sys.getenv("SIGNAL_AXIS_PREFIX", "signal_axis_map_phylo")   # e.g. ..._16S_fallback_test to browse the 16S tree
TIPS_F <- paste0(PREFIX, "_tips.csv")
MEM_F  <- "taxa_axis_membership.csv"
NWK_F  <- paste0(PREFIX, ".nwk")
PROV_F <- paste0(PREFIX, "_provenance.json")

find_file <- function(f, required = TRUE) {
  for (d in c("data", ".", "..", file.path("..", "data"))) {
    p <- file.path(d, f)
    if (file.exists(p)) return(p)
  }
  if (required) stop(sprintf("Cannot find %s. Put %s, %s and %s in a 'data' folder next to this script.",
                             f, TIPS_F, MEM_F, NWK_F), call. = FALSE)
  NA_character_
}

# ---- Newick (no ape): flat parent/length/label vectors; a parent's index < its children's
parse_newick <- function(s) {
  s <- sub(";$", "", gsub("[[:space:]]", "", paste(s, collapse = "")))
  tok <- regmatches(s, gregexpr("[(),:]|[^(),:]+", s))[[1]]
  parent <- integer(length(tok)); blen <- numeric(length(tok)); label <- rep(NA_character_, length(tok))
  n <- 0L; stack <- integer(0); last <- NA_integer_; want <- FALSE
  for (t in tok) {
    if (t == "(") {
      n <- n + 1L; parent[n] <- if (length(stack)) stack[length(stack)] else 0L
      stack <- c(stack, n); last <- NA_integer_
    } else if (t == ",") {
      last <- NA_integer_
    } else if (t == ")") {
      last <- stack[length(stack)]; stack <- stack[-length(stack)]
    } else if (t == ":") {
      want <- TRUE
    } else if (want) {
      blen[last] <- as.numeric(t); want <- FALSE
    } else if (is.na(last)) {
      n <- n + 1L; parent[n] <- stack[length(stack)]; label[n] <- t; last <- n
    } else {
      label[last] <- t
    }
  }
  if (length(stack)) stop("unbalanced Newick string")
  list(parent = parent[seq_len(n)], blen = blen[seq_len(n)], label = label[seq_len(n)])
}

# ---- circular layout: same rules as layout_tree() in microbe_phylogeny.py
#      (ladderised: smaller clades first, ties by C-locale label; tips evenly spaced;
#       internal node angle = midpoint of first/last child; radius = sqrt or linear depth)
layout_tree <- function(tr, radial = "sqrt", gap_deg = 6, R_TIP = 6) {
  n <- length(tr$parent); par <- tr$parent
  kids <- split(seq_len(n)[par > 0], par[par > 0])
  ch <- vector("list", n); ch[as.integer(names(kids))] <- kids
  ntips <- integer(n); minlab <- character(n); ordered <- vector("list", n)
  for (i in n:1) {
    k <- ch[[i]]
    if (is.null(k)) { ntips[i] <- 1L; minlab[i] <- tr$label[i]
    } else {
      ntips[i] <- sum(ntips[k]); minlab[i] <- sort(minlab[k], method = "radix")[1]
      ordered[[i]] <- k[order(ntips[k], minlab[k], method = "radix")]
    }
  }
  tip_order <- integer(0); stack <- 1L
  while (length(stack)) {
    i <- stack[length(stack)]; stack <- stack[-length(stack)]
    oc <- ordered[[i]]
    if (is.null(oc)) tip_order <- c(tip_order, i) else stack <- c(stack, rev(oc))
  }
  step <- (2 * pi - gap_deg * pi / 180) / length(tip_order)
  theta <- numeric(n)
  theta[tip_order] <- gap_deg * pi / 360 + (seq_along(tip_order) - 0.5) * step
  for (i in n:1) if (!is.null(ordered[[i]])) theta[i] <- (theta[ordered[[i]][1]] + theta[ordered[[i]][length(ordered[[i]])]]) / 2
  d <- numeric(n)
  for (i in seq_len(n)) d[i] <- (if (par[i] > 0) d[par[i]] else 0) + tr$blen[i]
  r <- if (radial == "sqrt") sqrt(d) else d
  r <- r / max(r[tip_order]) * R_TIP
  list(theta = theta, r = r, depth = d, tip_order = tip_order, step = step)
}

# ---- data
TIPS <- read_csv(find_file(TIPS_F), show_col_types = FALSE)
MEM  <- read_csv(find_file(MEM_F),  show_col_types = FALSE)
TR   <- parse_newick(readLines(find_file(NWK_F), warn = FALSE))
pf   <- find_file(PROV_F, required = FALSE)
PROV <- if (!is.na(pf)) paste(readLines(pf, warn = FALSE), collapse = " ") else ""
ROUTE <- if (grepl('"route": "gtdb"', PROV)) "gtdb" else if (grepl('"route": "16S"', PROV)) "16S" else "unknown"
REL   <- sub('.*"gtdb_version": "([^"]+)".*', "\\1", PROV); if (!grepl("^v", REL)) REL <- ""

for (cc in c("gtdb_species", "gtdb_lineage", "source_name", "purity", "n_genomes", "rule"))
  if (!cc %in% names(TIPS)) TIPS[[cc]] <- NA
for (cc in c("renamed", "gtdb_name_differs"))
  if (!cc %in% names(TIPS)) TIPS[[cc]] <- FALSE

N <- length(TR$parent)
TIP_NODES <- which(!(seq_len(N) %in% TR$parent))
TIP <- tibble(node = TIP_NODES, source_id = TR$label[TIP_NODES]) |>
  left_join(TIPS, by = "source_id")
if (anyNA(TIP$taxid)) stop("Tree tips missing from ", TIPS_F, ": ",
                           paste(TIP$source_id[is.na(TIP$taxid)], collapse = ", "))
TIP$phylum <- ifelse(is.na(TIP$gtdb_lineage), NA_character_, sub(".*p__([^;]*).*", "\\1", TIP$gtdb_lineage))
UNPLACED <- sum(is.na(TIPS$source_id) | !TIPS$source_id %in% TIP$source_id)

AX <- MEM |> distinct(axis_idx, axis, direction, axis_color) |> arrange(axis_idx)
MEMT <- MEM |> distinct(taxid, axis_idx)
TIP$axes_n   <- as.integer(table(factor(MEMT$taxid, levels = TIP$taxid))[as.character(TIP$taxid)])
TIP$axes_str <- vapply(TIP$taxid, function(t) paste(AX$axis[AX$axis_idx %in% MEMT$axis_idx[MEMT$taxid == t]], collapse = "; "), "")

abbr <- function(x) ifelse(is.na(x), "", sub("^(\\S)\\S*", "\\1.", sub("\\s*\\(.*$", "", x)))
TIP$txt <- {
  base <- ifelse(TIP$kind == "species",
                 paste0(TIP$label, ifelse(TIP$renamed %in% TRUE, "\u2020", ""), ifelse(TIP$gtdb_name_differs %in% TRUE, "\u2021", "")),
                 ifelse(TIP$rank == "genus", paste0(TIP$label, " spp."), paste0(TIP$label, " (", TIP$rank, ")")))
  repn <- abbr(dplyr::coalesce(TIP$gtdb_species, TIP$source_name))
  paste0(base, ifelse(TIP$kind == "clade_rep" & nzchar(repn), paste0(" [", repn, "]"), ""))
}
TIP$face <- ifelse(TIP$kind == "species" | TIP$rank == "genus", "italic", "plain")

# arbitrary stems (GTDB route): the two root branches join separately-inferred trees
ARB <- if (ROUTE == "gtdb") which(TR$parent == 1L) else integer(0)
ALL_EDGES <- which(TR$parent > 0)
TOPCHILD <- integer(N); for (i in 2:N) TOPCHILD[i] <- if (TR$parent[i] == 1L) i else TOPCHILD[TR$parent[i]]
GROUP <- if (ROUTE == "gtdb") TOPCHILD else rep(1L, N)      # which separately inferred tree a node belongs to
TOTAL_LEN <- sum(TR$blen[setdiff(ALL_EDGES, ARB)])

# ---- geometry (rings, cells, edges, tips) ------------------------------------------------
R_TIP <- 6; DR <- 0.15; GAP <- 0.2
ring_r <- numeric(nrow(AX)); cur <- R_TIP + 0.5
for (k in seq_len(nrow(AX))) {
  if (k > 1 && AX$direction[k] != AX$direction[k - 1]) cur <- cur + GAP
  ring_r[k] <- cur; cur <- cur + DR
}
R_OUT <- max(ring_r) + DR
KR <- nrow(AX); KT <- nrow(TIP)

edge_paths <- function(L) {
  nodes <- which(TR$parent > 0)
  bind_rows(lapply(nodes, function(i) {
    p <- TR$parent[i]; th <- seq(L$theta[p], L$theta[i], length.out = 14)
    tibble(node = i, x = c(L$r[p] * sin(th), L$r[i] * sin(L$theta[i])),
           y = c(L$r[p] * cos(th), L$r[i] * cos(L$theta[i])))
  }))
}
make_geom <- function(radial) {
  L <- layout_tree(TR, radial, R_TIP = R_TIP)
  th <- L$theta[TIP$node]; rr <- L$r[TIP$node]; deg <- (th * 180 / pi) %% 360
  tips <- TIP |> mutate(theta = th, x = rr * sin(th), y = rr * cos(th),
                        lx = (R_OUT + 0.12) * sin(th), ly = (R_OUT + 0.12) * cos(th),
                        angle = ifelse(deg <= 180, 90 - deg, 270 - deg), hj = ifelse(deg <= 180, 0, 1))
  guides <- tibble(x0 = tips$x, y0 = tips$y, x1 = (R_TIP + 0.45) * sin(th), y1 = (R_TIP + 0.45) * cos(th))
  list(L = L, edges = edge_paths(L), tips = tips, guides = guides)
}
GEOM <- list(sqrt = make_geom("sqrt"), linear = make_geom("linear"))

# ring cells: one 12-vertex annular sector per (tip, axis); group index gi = (tip-1)*KR + ring
CI <- expand.grid(ring = seq_len(KR), tip = seq_len(KT))
CI$axis_idx <- AX$axis_idx[CI$ring]; CI$color <- AX$axis_color[CI$ring]
CI$member <- paste(TIP$taxid[CI$tip], CI$axis_idx) %in% paste(MEMT$taxid, MEMT$axis_idx)
CELLP <- local({
  g <- expand.grid(v = 1:12, ring = seq_len(KR), tip = seq_len(KT))
  rel <- c(seq(-0.5, 0.5, length.out = 6), seq(0.5, -0.5, length.out = 6))
  th <- GEOM$sqrt$tips$theta[g$tip] + rel[g$v] * GEOM$sqrt$L$step * 0.92
  r  <- ring_r[g$ring] + ifelse(g$v <= 6, DR * 0.92, 0)
  data.frame(gi = (g$tip - 1) * KR + g$ring, x = r * sin(th), y = r * cos(th))
})

# ---- selection logic ---------------------------------------------------------------------
select_rows <- function(sel, match = "any") {
  if (!length(sel)) return(integer(0))
  m <- MEMT[MEMT$axis_idx %in% sel, ]
  k <- table(m$taxid)
  keep <- if (match == "all") as.numeric(names(k)[k == length(sel)]) else as.numeric(names(k))
  which(TIP$taxid %in% keep)
}
span_edges <- function(nodes) {          # edges of the smallest subtree connecting `nodes`
  m <- length(nodes)
  if (m == 0) return(integer(0))
  if (m == 1) return(nodes)              # single taxon: show its own branch (breadth stays 0)
  cnt <- integer(N); cnt[nodes] <- 1L
  for (i in N:2) cnt[TR$parent[i]] <- cnt[TR$parent[i]] + cnt[i]
  which(cnt >= 1L & cnt < m & TR$parent > 0)
}
highlight_edges <- function(nodes) {      # per separately-inferred tree, so the join is never highlighted
  if (!length(nodes)) return(integer(0))
  # a lone taxon in its own tree (the archaeon) has the arbitrary stem as its terminal branch: never draw it
  setdiff(unique(unlist(lapply(split(nodes, GROUP[nodes]), span_edges))), ARB)
}
breadth_of <- function(rows) {
  if (length(rows) < 2) return(0)
  nodes <- TIP$node[rows]
  sum(vapply(split(nodes, GROUP[nodes]), function(nd)
    if (length(nd) < 2) 0 else sum(TR$blen[setdiff(span_edges(nd), ARB)]), numeric(1)))
}
summarise_rows <- function(rows) {
  d <- TIP[rows, ]; b <- breadth_of(rows)
  list(n = nrow(d), n_sp = sum(d$kind == "species"), n_rep = sum(d$kind == "clade_rep"),
       phyla = sort(unique(stats::na.omit(d$phylum))), breadth = b, pct = 100 * b / TOTAL_LEN,
       cross = length(unique(GROUP[d$node])) > 1)
}
fmt_len <- function(x) formatC(x, digits = 3, format = "g")

OVERVIEW <- {
  o <- lapply(AX$axis_idx, function(a) {
    s <- summarise_rows(select_rows(a)); 
    data.frame(Axis = paste0(a, "  ", AX$axis[AX$axis_idx == a]),
               Direction = ifelse(AX$direction[AX$axis_idx == a] == "microbe_to_neuron", "microbe \u2192 neuron", "neuron \u2192 microbe"),
               Taxa = s$n, Species = s$n_sp, `Clade reps` = s$n_rep, Phyla = length(s$phyla),
               Breadth = s$breadth, `% of tree` = round(s$pct, 1), check.names = FALSE)
  })
  o <- bind_rows(o); o[order(-o$Breadth, o$Axis), ]
}
OVERVIEW$Breadth <- as.numeric(fmt_len(OVERVIEW$Breadth))

# ---- the plot ------------------------------------------------------------------------------
build_plot <- function(sel, rows, labels = "highlighted", radial = "sqrt", muted = TRUE,
                       branches = TRUE, clicked = NULL) {
  G <- GEOM[[radial]]; has <- length(sel) > 0
  col1 <- if (length(sel) == 1) AX$axis_color[AX$axis_idx == sel] else "#1B1B1B"
  hl_edges <- if (has && branches) highlight_edges(TIP$node[rows]) else integer(0)

  E <- G$edges; E$hl <- E$node %in% hl_edges; E$arb <- E$node %in% ARB
  E$col <- ifelse(E$hl, col1, ifelse(E$arb, "#C9959B", "grey68"))
  E$lw  <- ifelse(E$hl, 1.0, 0.35); E$lt <- ifelse(E$arb, "22", "solid")

  fill <- ifelse(CI$member, CI$color, "#F2F2F2"); alpha <- rep(1, nrow(CI))
  if (has) {
    insel <- CI$axis_idx %in% sel; tip_in <- (seq_len(KT) %in% rows)[CI$tip]
    alpha <- ifelse(!CI$member, 1, ifelse(insel & tip_in, 1, ifelse(insel, 0.35, ifelse(muted, 0.2, 0))))
    fill  <- ifelse(CI$member & (insel | muted), CI$color, "#F2F2F2")
    alpha[fill == "#F2F2F2"] <- 1
  }
  P <- CELLP; P$fill <- fill[P$gi]; P$alpha <- alpha[P$gi]

  T <- G$tips; T$hl <- seq_len(KT) %in% rows; is_rep <- T$kind == "clade_rep"
  if (!has) { T$f <- ifelse(is_rep, "white", "black"); T$o <- "black"; T$sz <- 2.2
  } else    { T$f <- ifelse(T$hl, ifelse(is_rep, "white", col1), "grey88"); T$o <- ifelse(T$hl, col1, "grey65")
              T$sz <- ifelse(T$hl, 3.6, 1.6) }
  T$lcol <- ifelse(!has | T$hl, "grey5", "grey60")
  Lb <- switch(labels, all = T, none = T[0, ], highlighted = if (has) T[T$hl, ] else T)

  lim <- if (nrow(Lb)) {
    rad <- R_OUT + 0.12 + 0.14 * nchar(Lb$txt)
    max(R_OUT + 0.35, abs(rad * sin(Lb$theta)) + 0.35, abs(rad * cos(Lb$theta)) + 0.35)
  } else R_OUT + 0.6
  g <- ggplot() +
    geom_polygon(data = P, aes(x, y, group = gi, fill = fill, alpha = alpha), colour = NA) +
    geom_segment(data = G$guides, aes(x = x0, y = y0, xend = x1, yend = y1), colour = "grey85",
                 linewidth = 0.25, linetype = "dotted") +
    geom_path(data = E[!E$hl, ], aes(x, y, group = node, colour = col, linewidth = lw, linetype = lt), lineend = "round") +
    geom_path(data = E[E$hl, ],  aes(x, y, group = node, colour = col, linewidth = lw, linetype = lt), lineend = "round") +
    geom_point(data = T, aes(x, y, fill = f, colour = o, size = sz), shape = 21, stroke = 0.6)
  if (nrow(Lb))
    g <- g + geom_text(data = Lb, aes(lx, ly, label = txt, angle = angle, hjust = hj, fontface = face, colour = lcol), size = 2.7)
  if (has) {
    ax_lab <- tibble(x = 0, y = ring_r[AX$axis_idx %in% sel] + DR / 2, lab = AX$axis_idx[AX$axis_idx %in% sel])
    g <- g + geom_text(data = ax_lab, aes(x, y, label = lab), size = 2.1, colour = "grey10", fontface = "bold")
  }
  if (!is.null(clicked) && clicked %in% T$source_id)
    g <- g + geom_point(data = T[T$source_id == clicked, ], aes(x, y), shape = 21, size = 7, fill = NA,
                        colour = "#E69F00", stroke = 1.4)
  g + scale_fill_identity() + scale_alpha_identity() + scale_colour_identity() + scale_linewidth_identity() +
    scale_linetype_identity() + scale_size_identity() +
    coord_fixed(xlim = c(-lim, lim), ylim = c(-lim, lim), clip = "off") +
    theme_void() + theme(plot.margin = margin(4, 4, 4, 4)) -> g
  attr(g, "lim") <- lim
  g
}

# ---- UI -------------------------------------------------------------------------------------
axis_choices <- function(dir) {
  a <- AX[AX$direction == dir, ]
  n_ax <- vapply(a$axis_idx, function(i) length(select_rows(i)), 1L)
  list(names = lapply(seq_len(nrow(a)), function(j)
         HTML(sprintf('<span class="chip" style="background:%s"></span>%d&nbsp; %s <span class="muted">(%d)</span>',
                      a$axis_color[j], a$axis_idx[j], a$axis[j], n_ax[j]))),
       values = as.character(a$axis_idx))
}
CH_M <- axis_choices("microbe_to_neuron"); CH_N <- axis_choices("neuron_to_microbe")
DEFAULT_AXIS <- as.character(AX$axis_idx[match("GABA", AX$axis)]); if (is.na(DEFAULT_AXIS)) DEFAULT_AXIS <- character(0)
src_line <- switch(ROUTE,
  gtdb = sprintf("GTDB %s reference trees, pruned. Branch length: substitutions/site (protein alignment).", REL),
  `16S` = "RefSeq 16S rRNA, MAFFT + FastTree. Branch length: substitutions/site (16S).",
  "Tree source unknown (provenance file not found).")

ui <- fluidPage(
  tags$head(tags$style(HTML("
    .chip{display:inline-block;width:11px;height:11px;border-radius:2px;margin-right:5px;vertical-align:-1px}
    .axbox{max-height:340px;overflow-y:auto;border:1px solid #ddd;border-radius:4px;padding:4px 10px;margin-bottom:10px}
    .axbox .checkbox{margin:2px 0}.muted{color:#888;font-size:85%}
    .metric{background:#f6f7f9;border-radius:4px;padding:8px 12px;margin-bottom:8px}
    .caveat{color:#8a6d3b;font-size:85%}"))),
  titlePanel("Signal axes on the microbial phylogeny"),
  sidebarLayout(
    sidebarPanel(width = 4,
      h5(strong("Signal axes"), span(class = "muted", "  (number of taxa)")),
      div(class = "axbox",
          checkboxGroupInput("ax_m2n", "Microbe \u2192 neuron", choiceNames = CH_M$names, choiceValues = CH_M$values,
                             selected = DEFAULT_AXIS),
          if (length(CH_N$values)) checkboxGroupInput("ax_n2m", "Neuron \u2192 microbe", choiceNames = CH_N$names,
                                                      choiceValues = CH_N$values)),
      radioButtons("match", "When several axes are selected, highlight taxa listed under",
                   c("any selected axis (union)" = "any", "all selected axes (intersection)" = "all")),
      radioButtons("labels", "Tip labels", c("highlighted only" = "highlighted", "all taxa" = "all", "none" = "none"),
                   inline = TRUE),
      radioButtons("radial", "Radius scale", c("square root of distance" = "sqrt", "linear (true scale)" = "linear"),
                   inline = TRUE),
      checkboxInput("branches", "Highlight connecting branches", TRUE),
      checkboxInput("muted", "Show other axes' rings, muted", TRUE),
      actionButton("clear", "Clear selection"),
      hr(), downloadButton("dl_csv", "Selected taxa (CSV)"), downloadButton("dl_png", "Tree (PNG)")),
    mainPanel(width = 8,
      uiOutput("summary"),
      tabsetPanel(
        tabPanel("Tree", plotOutput("tree", height = "820px", click = "tree_click"), uiOutput("tip_info")),
        tabPanel("Selected taxa", br(), tableOutput("tbl")),
        tabPanel("All axes", br(), p(class = "muted", "Sorted by breadth. Computed on the unrooted subtree connecting each axis's taxa."),
                 tableOutput("overview")),
        tabPanel("About", br(),
          p(src_line),
          p("Filled tip = species; hollow tip = a named genus, family, class or phylum drawn at a representative species ",
            "(bracketed in the label). Rings run inner to outer in curation order; the gap separates neuron \u2192 microbe axes."),
          p("\u2020 current NCBI name differs from the curated name. \u2021 GTDB name differs from NCBI."),
          p("Dashed pale-red branches are arbitrary stems that join GTDB's separately inferred bacterial and archaeal trees; ",
            "they are excluded from breadth."),
          p(class = "caveat", "These are representative taxa named in panel_curation.py, not taxa measured in any cohort. ",
            "A taxon appearing under an axis means the curation lists it as an example, not that it was detected."),
          if (UNPLACED > 0) p(class = "caveat", sprintf("%d curated taxa could not be placed in this tree and are not shown.", UNPLACED))))))
)

# ---- server --------------------------------------------------------------------------------
server <- function(input, output, session) {
  clicked <- reactiveVal(NULL)
  sel  <- reactive(as.integer(c(input$ax_m2n, input$ax_n2m)))
  rows <- reactive(select_rows(sel(), input$match))
  plot_obj <- reactive(build_plot(sel(), rows(), input$labels, input$radial, isTRUE(input$muted),
                                  isTRUE(input$branches), clicked()))

  output$tree <- renderPlot(plot_obj(), res = 96)

  observeEvent(input$tree_click, {
    hit <- nearPoints(GEOM[[input$radial]]$tips, input$tree_click, xvar = "x", yvar = "y", threshold = 14, maxpoints = 1)
    clicked(if (nrow(hit)) hit$source_id else NULL)
  })
  observeEvent(input$clear, {
    updateCheckboxGroupInput(session, "ax_m2n", selected = character(0))
    updateCheckboxGroupInput(session, "ax_n2m", selected = character(0))
    clicked(NULL)
  })

  output$summary <- renderUI({
    s <- sel()
    if (!length(s)) return(div(class = "metric", "Select one or more signal axes to highlight their taxa on the tree."))
    nm <- paste(AX$axis[AX$axis_idx %in% s], collapse = " + ")
    r <- rows()
    if (!length(r)) return(div(class = "metric", strong(nm), br(),
                               "No taxon is listed under all of the selected axes. Switch to union, or drop an axis."))
    m <- summarise_rows(r)
    div(class = "metric",
        strong(nm), if (length(s) > 1) span(class = "muted", if (input$match == "all") "  (intersection)" else "  (union)"), br(),
        sprintf("%d taxa (%d species, %d clade representatives) \u00b7 %d phyla%s", m$n, m$n_sp, m$n_rep, length(m$phyla),
                if (length(m$phyla)) paste0(": ", paste(m$phyla, collapse = ", ")) else ""), br(),
        if (m$n >= 2) sprintf("Breadth %s (%.1f%% of the tree's %s total)", fmt_len(m$breadth), m$pct, fmt_len(TOTAL_LEN))
        else "Breadth 0 (a single taxon)",
        if (m$cross) span(class = "caveat", br(), "Selection spans both separately inferred trees (bacterial and archaeal); ",
                          "each is scored within its own tree, so the archaeon adds 0 and the dashed join is not counted."))
  })

  sel_table <- reactive({
    r <- rows()
    if (!length(r)) return(NULL)
    d <- TIP[r, ]
    data.frame(Taxon = d$txt, Rank = d$rank, Kind = d$kind, `GTDB species` = d$gtdb_species, Phylum = d$phylum,
               `N axes` = d$axes_n, Axes = d$axes_str, check.names = FALSE)[order(d$phylum, d$label), ]
  })
  output$tbl <- renderTable({
    t <- sel_table(); if (is.null(t)) data.frame(Message = "No taxa selected.") else t
  }, striped = TRUE, spacing = "s", na = "")
  output$overview <- renderTable(OVERVIEW, striped = TRUE, spacing = "s", digits = 3)

  output$tip_info <- renderUI({
    id <- clicked(); if (is.null(id)) return(p(class = "muted", "Click a tip to see its details."))
    d <- TIP[TIP$source_id == id, ]
    chips <- lapply(which(AX$axis_idx %in% MEMT$axis_idx[MEMT$taxid == d$taxid]), function(j)
      HTML(sprintf('<span class="chip" style="background:%s"></span>%d %s&nbsp;&nbsp; ', AX$axis_color[j], AX$axis_idx[j], AX$axis[j])))
    div(class = "metric", strong(d$txt), br(),
        if (!is.na(d$gtdb_species)) span("GTDB species: ", d$gtdb_species, br()),
        if (!is.na(d$gtdb_lineage)) span(class = "muted", d$gtdb_lineage, br()),
        sprintf("%s \u00b7 chosen by: %s", d$kind, d$rule),
        if (!is.na(d$purity)) sprintf(" \u00b7 %s genomes, %.0f%% in this GTDB species", format(d$n_genomes, big.mark = ","), 100 * d$purity),
        br(), chips)
  })

  output$dl_csv <- downloadHandler(filename = function() "selected_taxa.csv",
                                   content = function(file) write_csv(if (is.null(sel_table())) data.frame() else sel_table(), file))
  output$dl_png <- downloadHandler(filename = function() "signal_axis_phylo.png",
                                   content = function(file) ggsave(file, plot_obj(), width = 11, height = 11, dpi = 200, bg = "white"))
}

shinyApp(ui, server)
