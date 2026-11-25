## CellChat v2 — UPS-high vs UPS-low malignant epithelium <-> TME
##
## DESIGN
##   One CellChat object; UPS-high and UPS-low are two cell groups next to
##   every TME cell type (single-dataset workflow). The merged multi-dataset
##   workflow is for separate biological conditions, so it is not used.

suppressPackageStartupMessages({
  library(CellChat)
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(ComplexHeatmap)
})

set.seed(42)
cat("CellChat version:", as.character(packageVersion("CellChat")), "\n")

source(here::here("config.R"))
result_path <- step_dir("11_CellChat")

## Figure export: every figure as PNG (300 dpi) + SVG (svglite, editable in Inkscape)
stopifnot("Install svglite for SVG output: install.packages('svglite')" =
            requireNamespace("svglite", quietly = TRUE))
save_gg <- function(p, name, width, height) {            # ggplot / patchwork objects
  ggsave(file.path(result_path, paste0(name, ".png")), p, width = width, height = height,
         dpi = 300, limitsize = FALSE)
  ggsave(file.path(result_path, paste0(name, ".svg")), p, width = width, height = height,
         device = svglite::svglite, limitsize = FALSE)
}
save_base <- function(draw_fun, name, width, height) {   # base graphics / ComplexHeatmap
  png(file.path(result_path, paste0(name, ".png")), width = width, height = height,
      units = "in", res = 300)
  draw_fun(); dev.off()
  svglite::svglite(file.path(result_path, paste0(name, ".svg")), width = width, height = height)
  draw_fun(); dev.off()
}

## 1. Load objects
##    seurat_harmony = full object (all cell types)
##    seu_mal        = CNV-confirmed malignant epithelium only
load(file.path(OBJ_DIR, "seurat_all_annotated.RData"))   # seurat_harmony (script 01)
load(file.path(OBJ_DIR, "seu_mal_UPS_groups.RData"))     # seu_mal (script 05)

# changing UPS-High to UPS-high and UPS-Low to UPS-low
seu_mal$UPS_group <- factor(
  dplyr::recode(as.character(seu_mal$UPS_group),
                "UPS-High" = "UPS-high",
                "UPS-Low"  = "UPS-low"),
  levels = c("UPS-high", "UPS-low")
)

table(seu_mal$UPS_group)



annot_col      <- "annotation"        # cell-type column in seurat_harmony
state_col      <- "UPS_group"         # AUCell UPS median-split column in seu_mal
patient_col    <- "orig.ident"        # patient / sample ID column in seurat_harmony
epi_annot_name <- "Epithelial cells"  # epithelial label in annot_col

stopifnot(annot_col   %in% colnames(seurat_harmony@meta.data),
          patient_col %in% colnames(seurat_harmony@meta.data),
          state_col   %in% colnames(seu_mal@meta.data))
stopifnot(epi_annot_name %in% unique(as.character(seurat_harmony@meta.data[[annot_col]])))

stopifnot(all(c("UPS-high", "UPS-low") %in% as.character(seu_mal@meta.data[[state_col]])))
print(table(seu_mal@meta.data[[state_col]]))

## 2. Label column on the full object
##    CNV-confirmed malignant cells -> "UPS-high" / "UPS-low"
##    all other cells               -> original annotation
##    epithelial cells not CNV-confirmed as malignant -> removed
state_labels <- setNames(as.character(seu_mal@meta.data[[state_col]]), colnames(seu_mal))
overlap_bc   <- intersect(colnames(seu_mal), colnames(seurat_harmony))
cat("\nMalignant cells matched into full object:", length(overlap_bc), "/", ncol(seu_mal), "\n")
if (length(overlap_bc) < ncol(seu_mal)) {
  warning("Some seu_mal barcodes are missing from seurat_harmony - check barcode prefixes.")
}

lab <- setNames(as.character(seurat_harmony@meta.data[[annot_col]]), colnames(seurat_harmony))
lab[overlap_bc] <- state_labels[overlap_bc]
seurat_harmony$cellchat_label <- lab

drop_epi <- lab == epi_annot_name
cat("Epithelial cells not CNV-confirmed as malignant (removed):", sum(drop_epi), "\n")
seu <- subset(seurat_harmony, cells = colnames(seurat_harmony)[!drop_epi])
# rm(seurat_harmony); gc()

## 3. Cell-group size filter = CellChat default (min.cells = 10)
min_cells  <- 10
lab_counts <- table(seu$cellchat_label)
print(lab_counts)
keep_labels <- names(lab_counts)[lab_counts >= min_cells]
dropped     <- setdiff(names(lab_counts), keep_labels)
if (length(dropped) > 0) {
  cat("Groups with <", min_cells, "cells (excluded):", paste(dropped, collapse = ", "), "\n")
}
stopifnot(all(c("UPS-high", "UPS-low") %in% keep_labels))
seu <- subset(seu, cells = colnames(seu)[seu$cellchat_label %in% keep_labels])

mal_states <- c("UPS-high", "UPS-low")
state_cols <- c("UPS-high" = "#B2182B", "UPS-low" = "#2166AC")   # red = UPS-high, blue = UPS-low
tme_labels <- sort(setdiff(keep_labels, mal_states))
cat("TME groups:", paste(tme_labels, collapse = ", "), "\n")

## 4. QC for reporting: group sizes and patient composition
mal_idx <- seu$cellchat_label %in% mal_states
cat("\nUPS-high n =", sum(seu$cellchat_label == "UPS-high"),
    " | UPS-low n =", sum(seu$cellchat_label == "UPS-low"), "\n")
comp <- table(patient = seu@meta.data[[patient_col]][mal_idx],
              state   = seu$cellchat_label[mal_idx])
write.csv(as.data.frame.matrix(comp), file.path(result_path, "UPS_group_by_patient_counts.csv"))
print(round(100 * prop.table(comp, margin = 2), 1))

## 5. CellChat inference (tutorial workflow, defaults)
seu <- JoinLayers(seu, assay = "RNA")
data.input <- LayerData(seu, assay = "RNA", layer = "data")

## CellChat needs log-normalised data, not raw counts (counts are integers)
vals <- if (inherits(data.input, "dgCMatrix")) head(data.input@x, 1e5) else as.numeric(data.input[1:200, 1:200])
stopifnot("RNA 'data' layer looks like raw counts - run NormalizeData first" = any(vals %% 1 != 0))

meta <- data.frame(
  labels  = factor(seu$cellchat_label, levels = c(mal_states, tme_labels)),
  samples = factor(seu@meta.data[[patient_col]]),
  row.names = colnames(seu)
)

cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")

CellChatDB <- CellChatDB.human
CellChatDB.use <- subsetDB(CellChatDB)   
cellchat@DB <- CellChatDB.use


cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, type = "triMean")
cellchat <- filterCommunication(cellchat, min.cells = min_cells)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

saveRDS(cellchat, file.path(result_path, "cellchat_UPS_singleobject.rds"))

## 6. Tables (subsetCommunication; significant only)
#  To make the tables and figures without rerunning CellChat, start here
source(here::here("config.R"))
result_path <- step_dir("11_CellChat")
cellchat    <- readRDS(file.path(result_path, "cellchat_UPS_singleobject.rds"))
mal_states  <- c("UPS-high", "UPS-low")
state_cols  <- c("UPS-high" = "#B2182B", "UPS-low" = "#2166AC")
tme_labels  <- sort(setdiff(levels(cellchat@idents), mal_states))

write.csv(subsetCommunication(cellchat),
          file.path(result_path, "LR_all_significant.csv"), row.names = FALSE)
write.csv(subsetCommunication(cellchat, sources.use = mal_states, targets.use = tme_labels),
          file.path(result_path, "LR_malignant_to_TME.csv"), row.names = FALSE)
write.csv(subsetCommunication(cellchat, sources.use = tme_labels, targets.use = mal_states),
          file.path(result_path, "LR_TME_to_malignant.csv"), row.names = FALSE)

## Pathway level (slot.name = "netP")
write.csv(subsetCommunication(cellchat, slot.name = "netP", sources.use = mal_states, targets.use = tme_labels),
          file.path(result_path, "pathway_malignant_to_TME.csv"), row.names = FALSE)
write.csv(subsetCommunication(cellchat, slot.name = "netP", sources.use = tme_labels, targets.use = mal_states),
          file.path(result_path, "pathway_TME_to_malignant.csv"), row.names = FALSE)

## Totals from the aggregated network (net$count = number of significant
## L-R pairs, net$weight = summed probability), malignant <-> TME edges
tot <- data.frame(
  state           = mal_states,
  outgoing_count  = rowSums(cellchat@net$count [mal_states, tme_labels, drop = FALSE]),
  outgoing_weight = rowSums(cellchat@net$weight[mal_states, tme_labels, drop = FALSE]),
  incoming_count  = colSums(cellchat@net$count [tme_labels, mal_states, drop = FALSE]),
  incoming_weight = colSums(cellchat@net$weight[tme_labels, mal_states, drop = FALSE])
)
write.csv(tot, file.path(result_path, "totals_malignant_TME.csv"), row.names = FALSE)
print(tot)

## 7. Aggregated network circle plots (tutorial)
groupSize <- as.numeric(table(cellchat@idents))

save_base(function() {
  par(mfrow = c(1, 2), xpd = TRUE)
  netVisual_circle(cellchat@net$count,  vertex.weight = groupSize, weight.scale = TRUE,
                   label.edge = FALSE, title.name = "Number of interactions")
  netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = TRUE,
                   label.edge = FALSE, title.name = "Interaction strength")
}, "circle_whole_network", width = 7, height = 4)

## Per-state circles, tutorial loop: edge.weight.max = max(mat)
plot_state_circles <- function(mat, measure_lab) {
  par(mfrow = c(2, 2), xpd = TRUE)
  for (s in mal_states) {   # signals sent by each state
    mat2 <- matrix(0, nrow(mat), ncol(mat), dimnames = dimnames(mat)); mat2[s, ] <- mat[s, ]
    netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = TRUE,
                     edge.weight.max = max(mat), title.name = paste0(s, " outgoing (", measure_lab, ")"))
  }
  for (s in mal_states) {   # signals received by each state
    mat2 <- matrix(0, nrow(mat), ncol(mat), dimnames = dimnames(mat)); mat2[, s] <- mat[, s]
    netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = TRUE,
                     edge.weight.max = max(mat), title.name = paste0(s, " incoming (", measure_lab, ")"))
  }
}
save_base(function() plot_state_circles(cellchat@net$count,  "number"),
          "circle_UPS_states_count",  width = 9.33, height = 9.33)
save_base(function() plot_state_circles(cellchat@net$weight, "strength"),
          "circle_UPS_states_weight", width = 9.33, height = 9.33)

## 8. Signalling roles (tutorial) + interaction totals per state
n_path <- length(cellchat@netP$pathways)
ht_h   <- max(8, 0.4 * n_path)   # cm, sized to number of pathways
ht_out <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", width = 8, height = 27)
ht_in  <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", width = 8, height = 27)
save_base(function() ComplexHeatmap::draw(ht_out + ht_in, ht_gap = grid::unit(0.5, "cm")),
          "signalingRole_heatmap_out_in", width = 11, height =13)

## Number of inferred interactions and interaction strength per state,
## Fig 4C layout (CellChat compareInteractions style): one bar per state,
## value printed on the bar. Values from net$count / net$weight for
## malignant <-> TME edges (totals_malignant_TME.csv). One row per
## direction, because outgoing and incoming can differ in opposite ways.
bar_panel <- function(dir_lab, meas_lab, values, digits = NULL) {
  d <- data.frame(state = factor(mal_states, levels = mal_states), value = values)
  d$label <- if (is.null(digits)) format(round(d$value)) else sprintf(paste0("%.", digits, "f"), d$value)
  ggplot(d, aes(x = state, y = value, fill = state)) +
    geom_col(width = 0.6) +
    geom_text(aes(label = label), vjust = -0.3, size = 3.2) +
    scale_fill_manual(values = state_cols) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    labs(x = NULL, y = meas_lab, title = dir_lab) +
    theme_classic(base_size = 10) +
    theme(legend.position = "none", plot.title = element_text(size = 10, hjust = 0.5))
}
## Number of inferred interactions only (strength stays in totals_malignant_TME.csv)

tot_long <- data.frame(
  state     = factor(rep(mal_states, 4), levels = mal_states),
  direction = factor(rep(c("Sent to TME", "Received from TME"), each = 2, times = 2),
                     levels = c("Sent to TME", "Received from TME")),
  measure   = factor(rep(c("Number of interactions", "Interaction strength"), each = 4),
                     levels = c("Number of interactions", "Interaction strength")),
  value     = c(tot$outgoing_count, tot$incoming_count, tot$outgoing_weight, tot$incoming_weight)
)
tot_long$label <- ifelse(tot_long$measure == "Number of interactions",
                         as.character(round(tot_long$value)), sprintf("%.2f", tot_long$value))

p_tot <- ggplot(tot_long, aes(x = direction, y = value, fill = state)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  geom_text(aes(label = label), position = position_dodge(width = 0.75), vjust = -0.4, size = 2.8) +
  facet_wrap(~ measure, scales = "free_y", strip.position = "left") +
  scale_fill_manual(values = state_cols) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = NULL, fill = NULL) +
  theme_classic(base_size = 9) +
  theme(strip.background = element_blank(), strip.placement = "outside",
        strip.text = element_text(size = 9), legend.position = "top")
save_gg(p_tot, "interactions_number_strength_UPS_states", width = 5, height = 2.5)

## Interaction heatmaps, Fig 4E layout: CellChat's own netVisual_heatmap on
## the single object. Rows = senders, columns = receivers; UPS-high and
## UPS-low are separate rows and columns next to every TME cell type.
## Top bars = total incoming per receiver, right bars = total outgoing per
## sender. Colour = value (number or strength), not a difference.
ht_cnt <- netVisual_heatmap(cellchat, measure = "count",  color.heatmap = "Reds")
ht_wt  <- netVisual_heatmap(cellchat, measure = "weight", color.heatmap = "Reds")
save_base(function() ComplexHeatmap::draw(ht_cnt + ht_wt, ht_gap = grid::unit(1, "cm")),
          "heatmap_interactions_count_weight", width = 8, height = 4)

## 9. L-R bubble plots (tutorial)
n_lr <- function(src, tgt) length(unique(subsetCommunication(cellchat, sources.use = src, targets.use = tgt)$interaction_name_2))
bub_w <- 3 + 0.3 * 2 * length(tme_labels)

gg_b_out <- netVisual_bubble(cellchat, sources.use = mal_states, targets.use = tme_labels, angle.x = 45, remove.isolate = TRUE)
save_gg(gg_b_out, "bubble_malignant_to_TME",
        width = 6, height =7)

gg_b_in <- netVisual_bubble(cellchat, sources.use = tme_labels, targets.use = mal_states, angle.x = 45, remove.isolate = TRUE)
save_gg(gg_b_in, "bubble_TME_to_malignant",
        width = 5.5, height = 10)


## 10. Pathway information flow table + selected-pathway figures
##  (a) Information flow = summed communication probability (CellChat's
##      definition) from netP$prob, same units for both states.
##      Table + relative/absolute flow figures (custom ggplot); descriptive.
##  (b) For the pathways discussed in the text: L-R bubble plots with
##      UPS-high and UPS-low side by side (netVisual_bubble, signaling =),
##      one per direction, plus optional violin plots of the pathway
##      genes per cell group (plotGeneExpression).
prob <- cellchat@netP$prob   # source x target x pathway (significant L-R only)
flow <- data.frame(
  pathway     = dimnames(prob)[[3]],
  High_to_TME = apply(prob["UPS-high", tme_labels, , drop = FALSE], 3, sum),
  Low_to_TME  = apply(prob["UPS-low",  tme_labels, , drop = FALSE], 3, sum),
  TME_to_High = apply(prob[tme_labels, "UPS-high", , drop = FALSE], 3, sum),
  TME_to_Low  = apply(prob[tme_labels, "UPS-low",  , drop = FALSE], 3, sum),
  row.names   = NULL
)
flow <- flow[rowSums(flow[, -1]) > 0, ]
flow <- flow[order(-(flow$High_to_TME + flow$Low_to_TME)), ]
write.csv(flow, file.path(result_path, "pathway_information_flow_UPS_states.csv"), row.names = FALSE)
print(flow)

## Information flow figure (layout of rankNet comparison, Fig 4F style):
## left = relative information flow (UPS-high share, stacked to 1),
## right = information flow (summed probability, linear axis).
## Same pathway order in both panels. No significance colouring: rankNet's
## paired test exists only for the multi-dataset workflow.
plot_flow_4F <- function(df, hi, lo, title) {
  d <- df[(df[[hi]] + df[[lo]]) > 0, c("pathway", hi, lo)]
  names(d) <- c("pathway", "UPS-high", "UPS-low")
  d$rel_high <- d$`UPS-high` / (d$`UPS-high` + d$`UPS-low`)
  lev <- d$pathway[order(d$rel_high, d$`UPS-high` + d$`UPS-low`)]   # most UPS-high-biased at top
  dl <- tidyr::pivot_longer(d[, c("pathway", "UPS-high", "UPS-low")], c("UPS-high", "UPS-low"),
                            names_to = "state", values_to = "flow")
  dl$pathway <- factor(dl$pathway, levels = lev)
  dl$state   <- factor(dl$state, levels = c("UPS-low", "UPS-high"))   # UPS-high starts at 0 / sits on top
  
  p_rel <- ggplot(dl, aes(x = flow, y = pathway, fill = state)) +
    geom_col(position = "fill", width = 0.75, orientation = "y") +
    geom_vline(xintercept = 0.5, linetype = 2, linewidth = 0.3) +
    scale_fill_manual(values = state_cols, breaks = mal_states) +
    scale_x_continuous(expand = c(0, 0)) +
    labs(x = "Relative information flow", y = NULL, fill = NULL, title = title) +
    theme_classic(base_size = 9)
  
  p_abs <- ggplot(dl, aes(x = flow, y = pathway, fill = state)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.75, orientation = "y") +
    scale_fill_manual(values = state_cols, breaks = mal_states) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(x = "Information flow", y = NULL, fill = NULL) +
    theme_classic(base_size = 9)
  
  (p_rel | p_abs) + plot_layout(guides = "collect")
}

n_out <- sum((flow$High_to_TME + flow$Low_to_TME) > 0)
n_in  <- sum((flow$TME_to_High + flow$TME_to_Low) > 0)
save_gg(plot_flow_4F(flow, "High_to_TME", "Low_to_TME", "Malignant -> TME"),
        "information_flow_malignant_to_TME", width = 5, height = 3)
save_gg(plot_flow_4F(flow, "TME_to_High", "TME_to_Low", "TME -> Malignant"),
        "information_flow_TME_to_malignant", width = 5, height =3.2)

## (b) Pathways for the focused bubble plots, chosen by a fixed rule (not by eye):
##     per direction, pathways are ranked by the absolute difference in
##     information flow between UPS-high and UPS-low, and the top n_show are
##     plotted. n_show is a display choice (state it in the legend); the full
##     ranking is saved as CSV and the full bubble plots (section 9) show all.
n_show <- 6
rank_diff <- function(hi, lo) {
  d <- data.frame(pathway = flow$pathway, UPS_High = flow[[hi]], UPS_Low = flow[[lo]])
  d <- d[(d$UPS_High + d$UPS_Low) > 0, ]
  d$abs_diff  <- abs(d$UPS_High - d$UPS_Low)
  d$higher_in <- ifelse(d$UPS_High > d$UPS_Low, "UPS-high",
                        ifelse(d$UPS_High < d$UPS_Low, "UPS-low", "equal"))
  d[order(-d$abs_diff), ]
}
rank_out <- rank_diff("High_to_TME", "Low_to_TME")
rank_in  <- rank_diff("TME_to_High", "TME_to_Low")
write.csv(rank_out, file.path(result_path, "pathway_difference_ranking_malignant_to_TME.csv"), row.names = FALSE)
write.csv(rank_in,  file.path(result_path, "pathway_difference_ranking_TME_to_malignant.csv"), row.names = FALSE)
print(head(rank_out, n_show)); print(head(rank_in, n_show))

pw_out <- head(rank_out$pathway, n_show)
pw_in  <- head(rank_in$pathway,  n_show)
pathways_show <- union(pw_out, pw_in)   # for the optional violin plots

## Focused bubble plots: UPS-high and UPS-low side by side, malignant <-> TME edges only
n_sel  <- function(src, tgt, pw) nrow(unique(subsetCommunication(cellchat, sources.use = src,
                                                                 targets.use = tgt, signaling = pw)[, c("ligand", "receptor")]))

if (length(pw_out) > 0) {
  gg_sel_out <- netVisual_bubble(cellchat, sources.use = mal_states, targets.use = tme_labels,
                                 signaling = pw_out, angle.x = 45 , remove.isolate = TRUE)
  save_gg(gg_sel_out, "bubble_selected_pathways_malignant_to_TME",
          width = 5.5, height = 5.5)
}
if (length(pw_in) > 0) {
  gg_sel_in <- netVisual_bubble(cellchat, sources.use = tme_labels, targets.use = mal_states,
                                signaling = pw_in, angle.x = 45 , remove.isolate = TRUE)
  save_gg(gg_sel_in, "bubble_selected_pathways_TME_to_malignant",
          width =5.5, height = 6)
}

