## 08  UCPS high- vs low-risk differential expression and GSEA
## Input : results/objects/187tcga_list_training_validation_Datasets.Rdata  (script 04c)
##         results/objects/tcga_risk_malignant_only_tcga187.RData          (script 06)
## Output: results/08_risk_group_DE_GSEA/  DEG tables, GSEA tables, figures
##
## Pipeline: expression matrix -> limma-trend DE -> volcano -> ranked list
##           -> GSEA over 5 MSigDB collections -> UPS check
##
suppressPackageStartupMessages({
  library(limma); library(clusterProfiler)
  library(msigdbr); library(ggplot2); library(ggrepel)
  library(dplyr); library(stringr)
})
set.seed(42)

## 0. CONSTANTS -- defined once, never redefined below
FDR_CUT  <- 0.05   ## Benjamini & Hochberg 1995
PADJ_CUT <- 0.05
LFC_CUT  <- 0.5    ## |log2FC| > 0.5. Affects ONLY the volcano colouring and
## DEG_significant.csv. GSEA uses the full ranked list and is
## not affected by it.

## GSEA gene-set size bounds. clusterProfiler::GSEA defaults, kept unchanged and
## stated here because they DO exclude sets: any pathway with < 10 or > 500
## measured genes is never tested. Report both numbers in Methods.
MIN_GS <- 10
MAX_GS <- 500

## display only -- these never touch a p-value
TOP_EACH   <- 15   ## dot panels: sets kept per direction
TOP_HALL   <- 15   ## bar panel: sets kept per direction
WRAP_CHARS <- 42   ## label wrap width, so long GOBP labels fit inside 180 mm
LAB_PT     <- 8
HDR_PT     <- 8    ## point size of the two direction headers
DOT_MIN    <- 1.6  ## smallest plotted point diameter, mm -- aesthetic only
DOT_MAX    <- 5.5  ## largest  plotted point diameter, mm -- aesthetic only

## Significance colour ramp for the GSEA panels (blue = weaker, red = stronger FDR)
GRAD_LO  <- "#08519C"   ## colour at the weakest significant FDR
GRAD_HI  <- "#B2182B"   ## colour at the strongest FDR

## figure sizes, inches, set by hand (7.09 in = 180 mm double-column width)
W_VOLC <- 4; H_VOLC <- 4
W_HALL <- 7.09; H_HALL <- 6.50
W_REAC <- 7.09; H_REAC <- 8.50
W_KEGG <- 7.09; H_KEGG <- 7.50
W_GOBP <- 7.09; H_GOBP <- 9.50
W_C6   <- 7.09; H_C6   <- 7.50

COL <- c("High-risk" = "#C0392B", "Low-risk" = "#1F77B4")

## UPS term-matching patterns. Two forms because MSigDB gs_name uses underscores
## while the plotting `label` column has them replaced by spaces.
UPS_PAT_SET <- "UBIQUIT|PROTEASOM|DEUBIQUIT|SUMOYL|NEDDYL|CULLIN|PROTEIN_CATABOLIC"
UPS_PAT_LAB <- "ubiquit|proteasom|deubiquit|sumoyl|neddyl|cullin|protein catabolic"

## The final signature and the SIGN of each gene's Cox coefficient. This is the
## only place the signature is written down.
UCPS_GENES <- c("RARRES1","ASPHD1","HELLS","EXOSC4","KRT16","PUF60")
UCPS_POS   <- c("ASPHD1","EXOSC4","KRT16","PUF60")   ## coef > 0, up in high-risk
UCPS_NEG   <- c("RARRES1","HELLS")                   ## coef < 0, up in low-risk
stopifnot(setequal(UCPS_GENES, c(UCPS_POS, UCPS_NEG)))

source(here::here("config.R"))
RData_path <- OBJ_DIR
plot_path  <- step_dir("08_risk_group_DE_GSEA")
file_path  <- plot_path
saved <- character(0)

## 1. EXPRESSION MATRIX
load(file.path(OBJ_DIR, "187tcga_list_training_validation_Datasets.Rdata"))
load(file.path(RData_path, "tcga_risk_malignant_only_tcga187.RData"))

tcga_expr    <- list_train_vali_Data[["tcga"]]
tcga_expr$ID <- substr(tcga_expr$ID, 1, 12)
tcga_expr    <- tcga_expr[!duplicated(tcga_expr$ID), ]
tcga_expr    <- merge(tcga_expr, tcga_risk[, c("ID","risk_group","RS")], by = "ID")
cat("Merged samples:", nrow(tcga_expr), "\n"); print(table(tcga_expr$risk_group))
stopifnot(nrow(tcga_expr) > 100)

## GROUP-LABEL GATE. factor(x, levels = c("Low","High")) returns all-NA if the
## column is coded any other way ("low"/"high", "Low-risk"/"High-risk", 0/1),
## and model.matrix() on an all-NA factor returns a ZERO-ROW design instead of
## erroring -- lmFit then fails with an unrelated dimension message. Checked
## here where the fix is obvious.
cat("risk_group levels found:", paste(sort(unique(as.character(tcga_expr$risk_group))),
                                      collapse = " / "), "\n")
stopifnot(setequal(as.character(tcga_expr$risk_group), c("Low","High")))

gene_cols <- setdiff(colnames(tcga_expr), c("ID","OS.time","OS","risk_group","RS"))
expr_mat  <- t(as.matrix(tcga_expr[, gene_cols]))
storage.mode(expr_mat) <- "double"
colnames(expr_mat) <- tcga_expr$ID

## Undo R's uniquify suffix (.1, .2) that data.frame() added to duplicated gene
## symbols. Reported, not silent: this regex would also truncate a genuine
## symbol ending in ".<digits>", so the count must look like a duplicate count,
## not like a large fraction of the matrix.
new_rn <- sub("\\.[0-9]+$", "", gene_cols)
cat(sprintf("uniquify suffixes stripped from %d of %d column names\n",
            sum(new_rn != gene_cols), length(gene_cols)))
rownames(expr_mat) <- new_rn

## Non-gene columns that slipped into gene_cols become all-NA rows here.
na_rows <- sum(rowSums(is.na(expr_mat)) == ncol(expr_mat))
cat(sprintf("all-NA rows (non-numeric columns mis-read as genes): %d\n", na_rows))
cat(sprintf("NA cells in matrix: %d of %d\n", sum(is.na(expr_mat)), length(expr_mat)))

rv       <- rowSums((expr_mat - rowMeans(expr_mat, na.rm = TRUE))^2, na.rm = TRUE)
ord      <- order(rv, decreasing = TRUE)
expr_mat <- expr_mat[ord, , drop = FALSE][!duplicated(rownames(expr_mat)[ord]), , drop = FALSE]
cat("unique symbols:", nrow(expr_mat), "\n")

cat(sprintf("matrix range: %.3f to %.3f\n",
            min(expr_mat, na.rm = TRUE), max(expr_mat, na.rm = TRUE)))
if (max(expr_mat, na.rm = TRUE) > 50)
  stop("matrix looks LINEAR-scale: log2(x+1) it first or logFC is uninterpretable")
if (min(expr_mat, na.rm = TRUE) < 0)
  stop("matrix has negative values: it is not log2(x+1), so the `> 0` detection rule below is invalid")


DET_FRAC <- min(table(tcga_expr$risk_group)) / nrow(tcga_expr)
rv       <- rowSums((expr_mat - rowMeans(expr_mat, na.rm = TRUE))^2, na.rm = TRUE)
expr_mat <- expr_mat[rowMeans(expr_mat > 0, na.rm = TRUE) >= DET_FRAC & rv > 0, , drop = FALSE]
cat(sprintf("smallest risk group: %d of %d samples -> detection threshold %.1f%%\n",
            min(table(tcga_expr$risk_group)), nrow(tcga_expr), 100 * DET_FRAC))
cat(sprintf("genes retained: %d\n", nrow(expr_mat)))
stopifnot(all(UCPS_GENES %in% rownames(expr_mat)))
save(expr_mat, file = file.path(RData_path, "expr_mat_tcga_filtered.RData"))

## 2. DIFFERENTIAL EXPRESSION  (limma-trend, moderated t)
## limma-trend rather than voom because the input is already log2 -- voom needs
## raw counts. trend = TRUE models the mean-variance relationship on the log
## scale; robust = TRUE down-weights hypervariable genes in the prior estimate
## (Phipson et al., Ann Appl Stat 2016). Both are the standard choices here.
stopifnot(identical(colnames(expr_mat), tcga_expr$ID))
group <- factor(tcga_expr$risk_group, levels = c("Low","High"))
stopifnot(!any(is.na(group)))
design <- model.matrix(~ group)
stopifnot(nrow(design) == ncol(expr_mat))
fit   <- eBayes(lmFit(expr_mat, design), trend = TRUE, robust = TRUE)
deg   <- topTable(fit, coef = 2, number = Inf, sort.by = "t")
deg$gene <- rownames(deg)

## `significance` is assigned before either CSV is written, so the full table
## records which genes were called.
deg$significance <- "NS"
deg$significance[deg$adj.P.Val < FDR_CUT & deg$logFC >  LFC_CUT] <- "Up"
deg$significance[deg$adj.P.Val < FDR_CUT & deg$logFC < -LFC_CUT] <- "Down"

deg_sig <- deg[deg$significance != "NS", ]
cat(sprintf("FDR < %.2f alone             : %d genes\n",
            FDR_CUT, sum(deg$adj.P.Val < FDR_CUT)))
cat(sprintf("FDR < %.2f and |log2FC| > %g : %d genes  [up in high %d / up in low %d]\n",
            FDR_CUT, LFC_CUT, nrow(deg_sig),
            sum(deg_sig$logFC > 0), sum(deg_sig$logFC < 0)))

lf <- setNames(deg$logFC, deg$gene)[UCPS_GENES]
print(data.frame(gene = UCPS_GENES, expected = ifelse(UCPS_GENES %in% UCPS_POS, "+", "-"),
                 logFC = round(unname(lf), 3)))
stopifnot(!any(is.na(lf)))                      ## all six measured
stopifnot(all(lf[UCPS_POS] > 0), all(lf[UCPS_NEG] < 0))
cat("DIRECTION GATE: PASS -- all 6 signature genes consistent with high-minus-low\n")

write.csv(deg,     file.path(file_path, "DEG_high_vs_low_full.csv"),        row.names = FALSE)
write.csv(deg_sig, file.path(file_path, "DEG_high_vs_low_significant.csv"), row.names = FALSE)
saved <- c(saved, file.path(file_path, c("DEG_high_vs_low_full.csv",
                                         "DEG_high_vs_low_significant.csv")))

## 3. VOLCANO
lab <- deg %>% filter(significance != "NS") %>% group_by(significance) %>%
  slice_max(order_by = abs(t), n = 15) %>% ungroup()

p_volc <- ggplot(deg, aes(logFC, -log10(adj.P.Val))) +
  geom_point(aes(colour = significance), alpha = 0.5, size = 1.1) +
  scale_colour_manual(values = c(Down = COL[["Low-risk"]], NS = "grey80",
                                 Up = COL[["High-risk"]]),
                      breaks = c("Up","Down","NS"),
                      labels = c("Up in high-risk","Up in low-risk","NS"), name = NULL) +
  geom_vline(xintercept = c(-LFC_CUT, LFC_CUT), linetype = "dashed", linewidth = 0.3) +
  geom_hline(yintercept = -log10(FDR_CUT),      linetype = "dashed", linewidth = 0.3) +
  geom_text_repel(data = lab, aes(label = gene), size = 3.1, box.padding = 0.4,
                  max.overlaps = Inf, segment.size = 0.3,
                  min.segment.length = 0, seed = 42) +
  labs(x = expression(log[2]~"fold change (high-risk / low-risk)"),
       y = expression(-log[10]~"(BH-adjusted "*italic(P)*")")) +
  theme_bw(base_size = 11) + 
  theme(
    legend.position = "top",
    panel.grid.major = element_blank(), # Removes major grid lines
    panel.grid.minor = element_blank()  # Removes minor grid lines
  )

ggsave(file.path(plot_path, "volcano_DEG.png"), p_volc,
       width = 4, height = 4, dpi = 300)
saved <- c(saved, file.path(plot_path, "volcano_DEG.png"))


## 4. RANKED LIST + GENE-SET RETRIEVAL + PRE-FLIGHT CHECKS
## Ranking metric is the moderated t statistic: it carries direction and
## magnitude, is variance-stabilised, and unlike logFC it does not put noisy
## low-expression genes at the extremes of the list.
ranked_genes <- sort(setNames(deg$t, deg$gene), decreasing = TRUE)

## fgsea requires unique, finite, named values
cat(sprintf("ranked list: n = %d | duplicated names = %d | non-finite = %d\n",
            length(ranked_genes), sum(duplicated(names(ranked_genes))),
            sum(!is.finite(ranked_genes))))
stopifnot(!any(duplicated(names(ranked_genes))), all(is.finite(ranked_genes)))

## msigdbr renamed category/subcategory -> collection/subcollection at v10, and
## renamed the returned column gs_subcat -> gs_subcollection at the same time.
MSIG_NEW <- "collection" %in% names(formals(msigdbr::msigdbr))
cat("msigdbr", as.character(packageVersion("msigdbr")),
    "-> using", if (MSIG_NEW) "collection/subcollection" else "category/subcategory", "\n")

## C2 fetched ONCE, then split into Reactome and KEGG
## The KEGG subcollection is "CP:KEGG_LEGACY" from msigdbr v10 / MSigDB 2023.2
## onward and "CP:KEGG" in earlier releases, so it is resolved from what the
## installed database contains. 
args_c2 <- list(species = "Homo sapiens")
args_c2[[if (MSIG_NEW) "collection" else "category"]] <- "C2"
c2 <- do.call(msigdbr::msigdbr, args_c2)
sc_col <- grep("^gs_subcat$|^gs_subcollection$", names(c2), value = TRUE)[1]
stopifnot(!is.na(sc_col))
sc_vals <- unique(c2[[sc_col]])

kegg_sub <- if ("CP:KEGG_LEGACY" %in% sc_vals) "CP:KEGG_LEGACY" else
  if ("CP:KEGG"        %in% sc_vals) "CP:KEGG"        else NA_character_
stopifnot(!is.na(kegg_sub), "CP:REACTOME" %in% sc_vals)
cat("KEGG subcollection resolved to:", kegg_sub,
    "| CP:KEGG* available:", paste(grep("^CP:KEGG", sc_vals, value = TRUE), collapse = ", "), "\n")

t2g_reac <- dplyr::select(c2[c2[[sc_col]] == "CP:REACTOME", ], gs_name, gene_symbol)
t2g_kegg <- dplyr::select(c2[c2[[sc_col]] == kegg_sub,      ], gs_name, gene_symbol)
stopifnot(nrow(t2g_reac) > 0, nrow(t2g_kegg) > 0)

## Hallmark, GOBP, C6
args_hall <- list(species = "Homo sapiens")
args_hall[[if (MSIG_NEW) "collection" else "category"]] <- "H"
t2g_hall <- dplyr::select(do.call(msigdbr::msigdbr, args_hall), gs_name, gene_symbol)

args_gobp <- list(species = "Homo sapiens")
args_gobp[[if (MSIG_NEW) "collection"    else "category"]]    <- "C5"
args_gobp[[if (MSIG_NEW) "subcollection" else "subcategory"]] <- "GO:BP"
t2g_gobp <- dplyr::select(do.call(msigdbr::msigdbr, args_gobp), gs_name, gene_symbol)

args_c6 <- list(species = "Homo sapiens")
args_c6[[if (MSIG_NEW) "collection" else "category"]] <- "C6"
t2g_c6 <- dplyr::select(do.call(msigdbr::msigdbr, args_c6), gs_name, gene_symbol)

stopifnot(nrow(t2g_hall) > 0, nrow(t2g_gobp) > 0, nrow(t2g_c6) > 0)
cat(sprintf("gene sets retrieved: Hallmark %d | Reactome %d | KEGG %d | GOBP %d | C6 %d\n",
            dplyr::n_distinct(t2g_hall$gs_name), dplyr::n_distinct(t2g_reac$gs_name),
            dplyr::n_distinct(t2g_kegg$gs_name), dplyr::n_distinct(t2g_gobp$gs_name),
            dplyr::n_distinct(t2g_c6$gs_name)))

## identifier space: msigdbr TERM2GENE is HGNC symbols. A mismatched ID space
## returns a valid-looking result object built on almost no overlap.
hall_sym <- unique(t2g_hall$gene_symbol)
ov <- length(intersect(names(ranked_genes), hall_sym))
cat(sprintf("symbol overlap with Hallmark: %d of %d Hallmark genes\n", ov, length(hall_sym)))
stopifnot(ov > 0.5 * length(hall_sym))   ## sanity floor, not a statistic

## 5. GSEA -- five MSigDB collections, each BH-corrected WITHIN itself
## 5a Hallmark -> bidirectional bar (the high-level high- vs low-risk contrast)
## 5b-5e       -> bidirectional dot plots, point size = leading-edge gene count
##
## pvalueCutoff = 1 so nothing is hidden; the FDR threshold is applied once, at
## plotting, via PADJ_CUT. eps = 0 stops fgsea clamping very small p-values at
## the machine-precision boundary. seed = TRUE makes the permutation
## reproducible alongside the set.seed(42) above.
##
## Direction is taken from sign(NES) only 

## 5a. Hallmark -- bidirectional bar
r_hall <- GSEA(ranked_genes, TERM2GENE = t2g_hall, pvalueCutoff = 1,
               minGSSize = MIN_GS, maxGSSize = MAX_GS, eps = 0, seed = TRUE, verbose = FALSE)
d_hall <- as.data.frame(r_hall@result)
d_hall$label      <- gsub("_", " ", sub("^HALLMARK_", "", d_hall$ID))
d_hall$collection <- "Hallmark"
write.csv(d_hall, file.path(file_path, "GSEA_Hallmark_full.csv"), row.names = FALSE)
cat(sprintf("Hallmark      tested %4d sets | %3d at FDR < %.2f (%3d high-risk, %3d low-risk)\n",
            nrow(d_hall), sum(d_hall$p.adjust < PADJ_CUT), PADJ_CUT,
            sum(d_hall$p.adjust < PADJ_CUT & d_hall$NES > 0),
            sum(d_hall$p.adjust < PADJ_CUT & d_hall$NES < 0)))

sig_hall <- d_hall[d_hall$p.adjust < PADJ_CUT, ]
if (nrow(sig_hall) == 0) {
  cat("Hallmark      no set at FDR <", PADJ_CUT, "-- no figure written\n")
} else {
  top_hall <- sig_hall %>%
    mutate(direction = ifelse(NES > 0, "High-risk", "Low-risk")) %>%
    group_by(direction) %>%
    slice_max(abs(NES), n = TOP_HALL, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(wrapped = str_wrap(label, WRAP_CHARS), negLogFDR = -log10(p.adjust))
  top_hall$row <- factor(top_hall$ID, levels = top_hall$ID[order(top_hall$NES)])
  
  ## headers are built from the data: a side with no bars gets no header
  hdr_hall <- na.omit(data.frame(
    x   = c(if (any(top_hall$NES < 0)) min(top_hall$NES) / 2 else NA_real_,
            if (any(top_hall$NES > 0)) max(top_hall$NES) / 2 else NA_real_),
    lab = c("Enriched in low-risk", "Enriched in high-risk"),
    stringsAsFactors = FALSE))
  
  cat(sprintf("Hallmark      drawn %2d of %3d significant (%2d high / %2d low) | FDR %.1e to %.1e | Inf: %d\n",
              nrow(top_hall), nrow(sig_hall), sum(top_hall$NES > 0), sum(top_hall$NES < 0),
              min(top_hall$p.adjust), max(top_hall$p.adjust),
              sum(is.infinite(top_hall$negLogFDR))))
  
  p_hall <- ggplot(top_hall, aes(x = NES, y = row, fill = negLogFDR)) +
    geom_col(width = 0.7) +
    geom_vline(xintercept = 0, linewidth = 0.35, colour = "grey25") +
    geom_text(data = hdr_hall, inherit.aes = FALSE,
              aes(x = x, y = nrow(top_hall) + 0.95, label = lab),
              size = HDR_PT / .pt, fontface = "bold") +
    scale_y_discrete(labels = setNames(top_hall$wrapped, top_hall$ID),
                     expand = expansion(add = c(0.6, 1.9))) +
    scale_fill_gradient(low = GRAD_LO, high = GRAD_HI,
                        name = expression(-log[10] * "(FDR)")) +
    labs(x = "Normalised enrichment score", y = NULL) +
    theme_bw(base_size = 9) +
    theme(axis.text.y        = element_text(size = LAB_PT),
          legend.position    = "right",
          legend.key.height  = unit(0.42, "cm"),
          legend.title       = element_text(size = 8))
  
  ggsave(file.path(plot_path, "GSEA_Hallmark_NES_bar.png"), p_hall,
         width = 6, height = 4, dpi = 300, limitsize = FALSE)
  
  ggsave(file.path(plot_path, "GSEA_Hallmark_NES_bar.svg"), p_hall,
         width = 6, height = 4, dpi = 300, limitsize = FALSE)
  
  
}

## 5b. Reactome -- bidirectional dot plot
r_reac <- GSEA(ranked_genes, TERM2GENE = t2g_reac, pvalueCutoff = 1,
               minGSSize = MIN_GS, maxGSSize = MAX_GS, eps = 0, seed = TRUE, verbose = FALSE)
d_reac <- as.data.frame(r_reac@result)
d_reac$label      <- gsub("_", " ", sub("^REACTOME_", "", d_reac$ID))
d_reac$collection <- "Reactome"
write.csv(d_reac, file.path(file_path, "GSEA_Reactome_full.csv"), row.names = FALSE)
cat(sprintf("Reactome      tested %4d sets | %3d at FDR < %.2f (%3d high-risk, %3d low-risk)\n",
            nrow(d_reac), sum(d_reac$p.adjust < PADJ_CUT), PADJ_CUT,
            sum(d_reac$p.adjust < PADJ_CUT & d_reac$NES > 0),
            sum(d_reac$p.adjust < PADJ_CUT & d_reac$NES < 0)))

sig_reac <- d_reac[d_reac$p.adjust < PADJ_CUT, ]
if (nrow(sig_reac) == 0) {
  cat("Reactome      no set at FDR <", PADJ_CUT, "-- no figure written\n")
} else {
  top_reac <- sig_reac %>%
    mutate(direction = ifelse(NES > 0, "High-risk", "Low-risk")) %>%
    group_by(direction) %>%
    slice_max(abs(NES), n = TOP_EACH, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(wrapped   = str_wrap(label, WRAP_CHARS),
           negLogFDR = -log10(p.adjust),
           le_count  = lengths(strsplit(core_enrichment, "/")))
  top_reac$row <- factor(top_reac$ID, levels = top_reac$ID[order(top_reac$NES)])
  
  hdr_reac <- na.omit(data.frame(
    x   = c(if (any(top_reac$NES < 0)) min(top_reac$NES) / 2 else NA_real_,
            if (any(top_reac$NES > 0)) max(top_reac$NES) / 2 else NA_real_),
    lab = c("Enriched in low-risk", "Enriched in high-risk"),
    stringsAsFactors = FALSE))
  
  cat(sprintf("Reactome      drawn %2d of %3d significant (%2d high / %2d low) | FDR %.1e to %.1e | LE %d-%d | Inf: %d\n",
              nrow(top_reac), nrow(sig_reac), sum(top_reac$NES > 0), sum(top_reac$NES < 0),
              min(top_reac$p.adjust), max(top_reac$p.adjust),
              min(top_reac$le_count), max(top_reac$le_count),
              sum(is.infinite(top_reac$negLogFDR))))
  
  p_reac <- ggplot(top_reac, aes(x = NES, y = row)) +
    geom_vline(xintercept = 0, linewidth = 0.35, colour = "grey25") +
    geom_point(aes(colour = negLogFDR, size = le_count)) +
    geom_text(data = hdr_reac, inherit.aes = FALSE,
              aes(x = x, y = nrow(top_reac) + 0.95, label = lab),
              size = HDR_PT / .pt, fontface = "bold") +
    scale_y_discrete(labels = setNames(top_reac$wrapped, top_reac$ID),
                     expand = expansion(add = c(0.7, 1.9))) +
    scale_colour_gradient(low = GRAD_LO, high = GRAD_HI,
                          name = expression(-log[10] * "(FDR)")) +
    scale_size_continuous(range = c(DOT_MIN, DOT_MAX), name = "Leading-edge\ngenes") +
    labs(x = "Normalised enrichment score", y = NULL) +
    theme_bw(base_size = 9) +
    theme(axis.text.y        = element_text(size = 8),
          panel.grid.major.y = element_line(linetype = "dotted", linewidth = 0.25,
                                            colour = "grey80"),
          legend.position    = "right",
          legend.key.height  = unit(0.42, "cm"),
          legend.title       = element_text(size = 8))
  ggsave(file.path(plot_path, "GSEA_Reactome_NES_dot.png"), p_reac,
         width = 7, height = 7, dpi = 300, limitsize = FALSE)
  ggsave(file.path(plot_path, "GSEA_Reactome_NES_dot.svg"), p_reac,
         width = W_REAC, height = H_REAC, dpi = 300, limitsize = FALSE)
}

## 5c. KEGG -- bidirectional dot plot
r_kegg <- GSEA(ranked_genes, TERM2GENE = t2g_kegg, pvalueCutoff = 1,
               minGSSize = MIN_GS, maxGSSize = MAX_GS, eps = 0, seed = TRUE, verbose = FALSE)
d_kegg <- as.data.frame(r_kegg@result)
d_kegg$label      <- gsub("_", " ", sub("^KEGG_", "", d_kegg$ID))
d_kegg$collection <- sub("^CP:", "", kegg_sub)   ## records which KEGG release was used
write.csv(d_kegg, file.path(file_path, "GSEA_KEGG_legacy_full.csv"), row.names = FALSE)
cat(sprintf("KEGG          tested %4d sets | %3d at FDR < %.2f (%3d high-risk, %3d low-risk)\n",
            nrow(d_kegg), sum(d_kegg$p.adjust < PADJ_CUT), PADJ_CUT,
            sum(d_kegg$p.adjust < PADJ_CUT & d_kegg$NES > 0),
            sum(d_kegg$p.adjust < PADJ_CUT & d_kegg$NES < 0)))

sig_kegg <- d_kegg[d_kegg$p.adjust < PADJ_CUT, ]
if (nrow(sig_kegg) == 0) {
  cat("KEGG          no set at FDR <", PADJ_CUT, "-- no figure written\n")
} else {
  top_kegg <- sig_kegg %>%
    mutate(direction = ifelse(NES > 0, "High-risk", "Low-risk")) %>%
    group_by(direction) %>%
    slice_max(abs(NES), n = TOP_EACH, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(wrapped   = str_wrap(label, WRAP_CHARS),
           negLogFDR = -log10(p.adjust),
           le_count  = lengths(strsplit(core_enrichment, "/")))
  top_kegg$row <- factor(top_kegg$ID, levels = top_kegg$ID[order(top_kegg$NES)])
  
  hdr_kegg <- na.omit(data.frame(
    x   = c(if (any(top_kegg$NES < 0)) min(top_kegg$NES) / 2 else NA_real_,
            if (any(top_kegg$NES > 0)) max(top_kegg$NES) / 2 else NA_real_),
    lab = c("Enriched\nin low-risk", "Enriched\nin high-risk"),
    stringsAsFactors = FALSE))
  
  cat(sprintf("KEGG          drawn %2d of %3d significant (%2d high / %2d low) | FDR %.1e to %.1e | LE %d-%d | Inf: %d\n",
              nrow(top_kegg), nrow(sig_kegg), sum(top_kegg$NES > 0), sum(top_kegg$NES < 0),
              min(top_kegg$p.adjust), max(top_kegg$p.adjust),
              min(top_kegg$le_count), max(top_kegg$le_count),
              sum(is.infinite(top_kegg$negLogFDR))))
  
  p_kegg <- ggplot(top_kegg, aes(x = NES, y = row)) +
    geom_vline(xintercept = 0, linewidth = 0.35, colour = "grey25") +
    geom_point(aes(colour = negLogFDR, size = le_count)) +
    geom_text(data = hdr_kegg, inherit.aes = FALSE,
              aes(x = x, y = nrow(top_kegg) + 0.95, label = lab),
              size = HDR_PT / .pt, fontface = "bold") +
    scale_y_discrete(labels = setNames(top_kegg$wrapped, top_kegg$ID),
                     expand = expansion(add = c(0.7, 1.9))) +
    scale_colour_gradient(low = GRAD_LO, high = GRAD_HI,
                          name = expression(-log[10] * "(FDR)")) +
    scale_size_continuous(range = c(DOT_MIN, DOT_MAX), name = "Leading-edge\ngenes") +
    labs(x = "Normalised enrichment score", y = NULL) +
    theme_bw(base_size = 9) +
    theme(axis.text.y        = element_text(size = LAB_PT),
          panel.grid.major.y = element_line(linetype = "dotted", linewidth = 0.25,
                                            colour = "grey80"),
          legend.position    = "right",
          legend.key.height  = unit(0.42, "cm"),
          legend.title       = element_text(size = 8))
  ggsave(file.path(plot_path, "GSEA_KEGG_legacy_NES_dot.png"), p_kegg,
         width = 6, height = 4, dpi = 300, limitsize = FALSE)
  ggsave(file.path(plot_path, "GSEA_KEGG_legacy_NES_dot.svg"), p_kegg,
         width = 6, height = 4, dpi = 300, limitsize = FALSE)
}

## 5d. GO Biological Process -- bidirectional dot plot
r_gobp <- GSEA(ranked_genes, TERM2GENE = t2g_gobp, pvalueCutoff = 1,
               minGSSize = MIN_GS, maxGSSize = MAX_GS, eps = 0, seed = TRUE, verbose = FALSE)
d_gobp <- as.data.frame(r_gobp@result)
d_gobp$label      <- gsub("_", " ", sub("^GOBP_", "", d_gobp$ID))
d_gobp$collection <- "GOBP"
write.csv(d_gobp, file.path(file_path, "GSEA_GOBP_full.csv"), row.names = FALSE)
cat(sprintf("GOBP          tested %4d sets | %3d at FDR < %.2f (%3d high-risk, %3d low-risk)\n",
            nrow(d_gobp), sum(d_gobp$p.adjust < PADJ_CUT), PADJ_CUT,
            sum(d_gobp$p.adjust < PADJ_CUT & d_gobp$NES > 0),
            sum(d_gobp$p.adjust < PADJ_CUT & d_gobp$NES < 0)))

sig_gobp <- d_gobp[d_gobp$p.adjust < PADJ_CUT, ]
if (nrow(sig_gobp) == 0) {
  cat("GOBP          no set at FDR <", PADJ_CUT, "-- no figure written\n")
} else {
  top_gobp <- sig_gobp %>%
    mutate(direction = ifelse(NES > 0, "High-risk", "Low-risk")) %>%
    group_by(direction) %>%
    slice_max(abs(NES), n = TOP_EACH, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(wrapped   = str_wrap(label, WRAP_CHARS),
           negLogFDR = -log10(p.adjust),
           le_count  = lengths(strsplit(core_enrichment, "/")))
  top_gobp$row <- factor(top_gobp$ID, levels = top_gobp$ID[order(top_gobp$NES)])
  
  hdr_gobp <- na.omit(data.frame(
    x   = c(if (any(top_gobp$NES < 0)) min(top_gobp$NES) / 2 else NA_real_,
            if (any(top_gobp$NES > 0)) max(top_gobp$NES) / 2 else NA_real_),
    lab = c("Enriched in low-risk", "Enriched in high-risk"),
    stringsAsFactors = FALSE))
  
  cat(sprintf("GOBP          drawn %2d of %3d significant (%2d high / %2d low) | FDR %.1e to %.1e | LE %d-%d | Inf: %d\n",
              nrow(top_gobp), nrow(sig_gobp), sum(top_gobp$NES > 0), sum(top_gobp$NES < 0),
              min(top_gobp$p.adjust), max(top_gobp$p.adjust),
              min(top_gobp$le_count), max(top_gobp$le_count),
              sum(is.infinite(top_gobp$negLogFDR))))
  
  p_gobp <- ggplot(top_gobp, aes(x = NES, y = row)) +
    geom_vline(xintercept = 0, linewidth = 0.35, colour = "grey25") +
    geom_point(aes(colour = negLogFDR, size = le_count)) +
    geom_text(data = hdr_gobp, inherit.aes = FALSE,
              aes(x = x, y = nrow(top_gobp) + 0.95, label = lab),
              size = HDR_PT / .pt, fontface = "bold") +
    scale_y_discrete(labels = setNames(top_gobp$wrapped, top_gobp$ID),
                     expand = expansion(add = c(0.7, 1.9))) +
    scale_colour_gradient(low = GRAD_LO, high = GRAD_HI,
                          name = expression(-log[10] * "(FDR)")) +
    scale_size_continuous(range = c(DOT_MIN, DOT_MAX), name = "Leading-edge\ngenes") +
    labs(x = "Normalised enrichment score", y = NULL) +
    theme_bw(base_size = 9) +
    theme(axis.text.y        = element_text(size = LAB_PT),
          panel.grid.major.y = element_line(linetype = "dotted", linewidth = 0.25,
                                            colour = "grey80"),
          legend.position    = "right",
          legend.key.height  = unit(0.42, "cm"),
          legend.title       = element_text(size = 8))
  ggsave(file.path(plot_path, "GSEA_GOBP_NES_dot.png"), p_gobp,
         width = 7, height = 9, dpi = 300, limitsize = FALSE)
  ggsave(file.path(plot_path, "GSEA_GOBP_NES_dot.svg"), p_gobp,
         width = W_GOBP, height = H_GOBP, dpi = 300, limitsize = FALSE)
}

## 5e. C6 oncogenic signatures -- bidirectional dot plot
r_c6 <- GSEA(ranked_genes, TERM2GENE = t2g_c6, pvalueCutoff = 1,
             minGSSize = MIN_GS, maxGSSize = MAX_GS, eps = 0, seed = TRUE, verbose = FALSE)
d_c6 <- as.data.frame(r_c6@result)
d_c6$label      <- gsub("_", " ", d_c6$ID)   ## C6 IDs carry no common prefix
d_c6$collection <- "C6_oncogenic"
write.csv(d_c6, file.path(file_path, "GSEA_C6_oncogenic_full.csv"), row.names = FALSE)
cat(sprintf("C6_oncogenic  tested %4d sets | %3d at FDR < %.2f (%3d high-risk, %3d low-risk)\n",
            nrow(d_c6), sum(d_c6$p.adjust < PADJ_CUT), PADJ_CUT,
            sum(d_c6$p.adjust < PADJ_CUT & d_c6$NES > 0),
            sum(d_c6$p.adjust < PADJ_CUT & d_c6$NES < 0)))

sig_c6 <- d_c6[d_c6$p.adjust < PADJ_CUT, ]
if (nrow(sig_c6) == 0) {
  cat("C6_oncogenic  no set at FDR <", PADJ_CUT, "-- no figure written\n")
} else {
  top_c6 <- sig_c6 %>%
    mutate(direction = ifelse(NES > 0, "High-risk", "Low-risk")) %>%
    group_by(direction) %>%
    slice_max(abs(NES), n = TOP_EACH, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(wrapped   = str_wrap(label, WRAP_CHARS),
           negLogFDR = -log10(p.adjust),
           le_count  = lengths(strsplit(core_enrichment, "/")))
  top_c6$row <- factor(top_c6$ID, levels = top_c6$ID[order(top_c6$NES)])
  
  hdr_c6 <- na.omit(data.frame(
    x   = c(if (any(top_c6$NES < 0)) min(top_c6$NES) / 2 else NA_real_,
            if (any(top_c6$NES > 0)) max(top_c6$NES) / 2 else NA_real_),
    lab = c("Enriched in low-risk", "Enriched in high-risk"),
    stringsAsFactors = FALSE))
  
  cat(sprintf("C6_oncogenic  drawn %2d of %3d significant (%2d high / %2d low) | FDR %.1e to %.1e | LE %d-%d | Inf: %d\n",
              nrow(top_c6), nrow(sig_c6), sum(top_c6$NES > 0), sum(top_c6$NES < 0),
              min(top_c6$p.adjust), max(top_c6$p.adjust),
              min(top_c6$le_count), max(top_c6$le_count),
              sum(is.infinite(top_c6$negLogFDR))))
  
  p_c6 <- ggplot(top_c6, aes(x = NES, y = row)) +
    geom_vline(xintercept = 0, linewidth = 0.35, colour = "grey25") +
    geom_point(aes(colour = negLogFDR, size = le_count)) +
    geom_text(data = hdr_c6, inherit.aes = FALSE,
              aes(x = x, y = nrow(top_c6) + 0.95, label = lab),
              size = HDR_PT / .pt, fontface = "bold") +
    scale_y_discrete(labels = setNames(top_c6$wrapped, top_c6$ID),
                     expand = expansion(add = c(0.7, 1.9))) +
    scale_colour_gradient(low = GRAD_LO, high = GRAD_HI,
                          name = expression(-log[10] * "(FDR)")) +
    scale_size_continuous(range = c(DOT_MIN, DOT_MAX), name = "Leading-edge\ngenes") +
    labs(x = "Normalised enrichment score", y = NULL) +
    theme_bw(base_size = 9) +
    theme(axis.text.y        = element_text(size = LAB_PT),
          panel.grid.major.y = element_line(linetype = "dotted", linewidth = 0.25,
                                            colour = "grey80"),
          legend.position    = "right",
          legend.key.height  = unit(0.42, "cm"),
          legend.title       = element_text(size = 8))
  ggsave(file.path(plot_path, "GSEA_C6_oncogenic_NES_dot.png"), p_c6,
         width = W_C6, height = H_C6, dpi = 300, limitsize = FALSE)
  ggsave(file.path(plot_path, "GSEA_C6_oncogenic_NES_dot.svg"), p_c6,
         width = W_C6, height = H_C6, dpi = 300, limitsize = FALSE)
  
}

## 5f. manifest -- only files that were actually written are recorded
saved <- c(saved,
           file.path(file_path, c("GSEA_Hallmark_full.csv", "GSEA_Reactome_full.csv",
                                  "GSEA_KEGG_legacy_full.csv", "GSEA_GOBP_full.csv",
                                  "GSEA_C6_oncogenic_full.csv")))
gsea_pngs <- file.path(plot_path, c("GSEA_Hallmark_NES_bar.png", "GSEA_Reactome_NES_dot.png",
                                    "GSEA_KEGG_legacy_NES_dot.png", "GSEA_GOBP_NES_dot.png",
                                    "GSEA_C6_oncogenic_NES_dot.png"))
saved <- c(saved, gsea_pngs[file.exists(gsea_pngs)])

## 6. UPS CHECK
## The signature genes come from UPS-High vs UPS-Low malignant cells. Whether the
## bulk risk score tracks UPS activity is therefore an open question, not an
## assumption, and it is measured here rather than inferred from the pathway
## figures.
stopifnot(identical(colnames(expr_mat), tcga_expr$ID))   ## re-asserted: far from section 2

## 6a. Are any of the six signature genes themselves UPS genes?
## Checked against the same gene sets the GSEA used.
ups_set_names <- unique(c(
  t2g_kegg$gs_name[grepl(UPS_PAT_SET, t2g_kegg$gs_name, ignore.case = TRUE)],
  t2g_reac$gs_name[grepl(UPS_PAT_SET, t2g_reac$gs_name, ignore.case = TRUE)],
  t2g_gobp$gs_name[grepl(UPS_PAT_SET, t2g_gobp$gs_name, ignore.case = TRUE)]))
ups_ann_genes <- unique(c(
  t2g_kegg$gene_symbol[t2g_kegg$gs_name %in% ups_set_names],
  t2g_reac$gene_symbol[t2g_reac$gs_name %in% ups_set_names],
  t2g_gobp$gene_symbol[t2g_gobp$gs_name %in% ups_set_names]))
ucps_in_ups <- intersect(UCPS_GENES, ups_ann_genes)
cat(sprintf("\nUPS-annotated gene sets matched: %d (%d unique genes)\n",
            length(ups_set_names), length(ups_ann_genes)))
cat("UCPS genes annotated to a UPS set:",
    if (length(ucps_in_ups)) paste(ucps_in_ups, collapse = ", ") else "none of the 6", "\n")

## 6b. Direction of the UPS-related sets that GSEA actually tested.
## NOTE for Methods: p.adjust was computed WITHIN each collection, so the counts
## below pool five separate BH families.
ups_all <- bind_rows(d_hall, d_reac, d_kegg, d_gobp, d_c6)
ups_tab <- ups_all[grepl(UPS_PAT_LAB, ups_all$label, ignore.case = TRUE),
                   c("collection","label","setSize","NES","pvalue","p.adjust")]
write.csv(ups_tab, file.path(file_path, "UPS_sets_in_GSEA.csv"), row.names = FALSE)
cat(sprintf("UPS-related sets tested: %d | significant: %d | of those, NES > 0: %d\n",
            nrow(ups_tab), sum(ups_tab$p.adjust < PADJ_CUT),
            sum(ups_tab$p.adjust < PADJ_CUT & ups_tab$NES > 0)))
print(head(ups_tab[order(ups_tab$p.adjust), ], 15))
saved <- c(saved, file.path(file_path, "UPS_sets_in_GSEA.csv"))

## 7. MANIFEST -- a silent ggsave failure must not pass unnoticed
cat("\n=== outputs ===\n")
print(data.frame(file = basename(saved), exists = file.exists(saved),
                 kb = round(file.size(saved) / 1024, 1)))
stopifnot(all(file.exists(saved)))
