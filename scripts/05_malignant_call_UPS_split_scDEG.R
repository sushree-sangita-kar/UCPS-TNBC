## 05  Malignant cell calls from inferCNV, UPS-High / UPS-Low median split
##     within malignant cells, single-cell DEGs, intersection with bulk DEGs
##
## Input : resources/malignant_classification_anchored.csv        (stored calls, see script 03)
##         results/objects/epithelial_cells.RData                  (script 01)
##         resources/UPS_AUCell_scores_epithelial_cells.csv        (stored scores, see script 02)
##         results/objects/bulk_degs_187.RData                     (script 04d)
## Output: results/objects/seu_mal_UPS_groups.RData   (object: seu_mal)
##         results/05_malignant_UPS_scDEG/  DEG tables, candidate genes, figures

source(here::here("config.R"))
result_path <- step_dir("05_malignant_UPS_scDEG")

## 1. Malignant calls: the stored calls used in the study (made in script 03)
res <- read.csv(file.path(RES_DIR, "malignant_classification_anchored.csv"))
cat("epithelial cells with a CNV call:", nrow(res), "| malignant:", sum(res$malignant), "\n")

library(Seurat)

## 2. Sanity check: res must have all 6 samples
stopifnot(exists("res"))
print(table(res$sample, res$malignant))          # 6 rows expected
res$malignant_lab <- ifelse(res$malignant, "Malignant", "Normal-like")

## 3. Epithelial object (script 01) + stored UPS AUCell scores + CNV labels
## The UPS_AUCell values are the exact per-cell scores used in the study, read
## from resources/ (see script 02 for how they were computed). 
load(file.path(OBJ_DIR, "epithelial_cells.RData"))   # epithelial_cells

auc <- read.csv(file.path(RES_DIR, "UPS_AUCell_scores_epithelial_cells.csv"))
stopifnot(setequal(auc$cell, colnames(epithelial_cells)))
epithelial_cells$UPS_AUCell <- auc$UPS_AUCell[match(colnames(epithelial_cells), auc$cell)]

epi <- epithelial_cells
epi$malignant <- res$malignant_lab[match(colnames(epi), res$cell)]
print(table(epi$malignant, useNA="ifany"))
print(table(epi$orig.ident, epi$malignant, useNA="ifany"))

# cells without a CNV call are not plotted
cells_called <- colnames(epi)[!is.na(epi$malignant)]   # 9,307 cells
p1 <- DimPlot(epi, group.by = "malignant", cells = cells_called)
ggsave(file.path(result_path,"malignant_dimplot.svg"), plot = p1, width = 5, height = 4, dpi = 300)

## 4. Malignant-only object; UPS median split within malignant cells

seu_mal <- subset(epi, subset = malignant == "Malignant")
cat("malignant-only cells:", ncol(seu_mal), "\n")     # expect ~7029


ups_col <- grep("UPS_AUCell", colnames(seu_mal@meta.data), value=TRUE)
ups_col <- setdiff(ups_col, "UPS_group")[1]           # the numeric score, not the label
cat("UPS score column:", ups_col, "\n")
stopifnot(!is.na(ups_col))
m <- seu_mal@meta.data[[ups_col]]
seu_mal$UPS_group <- ifelse(m > median(m), "UPS-High", "UPS-Low")
cat("\n--- malignant-only UPS split by patient (must stay multi-patient) ---\n")
print(table(seu_mal$UPS_group, seu_mal$orig.ident))

p2<- FeaturePlot(epi, features = "UPS_AUCell") +
  scale_colour_viridis_c(name = "AUCell Score", option = "viridis") +
  ggtitle("Cell scoring") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))


ggsave(file.path(result_path,"ups_score_featureplot_all_epi.png"), plot = p2, width = 4, height = 3, dpi = 300)



p3<- DimPlot(seu_mal, group.by = "UPS_group")
ggsave(file.path(result_path,"UPS_group_dimplot.svg"), plot = p3, width = 4, height = 3, dpi = 300)


## 5. UPS-High vs UPS-Low DEGs in malignant cells (Wilcoxon; Seurat p_val_adj)
Idents(seu_mal) <- "UPS_group"
deg_mal <- FindMarkers(seu_mal, ident.1="UPS-High", ident.2="UPS-Low",
                       test.use="wilcox", logfc.threshold=0.25, min.pct=0.1)
deg_mal$gene <- rownames(deg_mal)
deg_sig <- deg_mal[deg_mal$p_val_adj < 0.05, ]
cat("\nmalignant-only DEG (padj<0.05):", nrow(deg_sig), "\n")

write.csv(deg_sig, file.path(result_path,"deg_sig_sc.csv"))

## 6. Intersect with TCGA bulk DEGs -> candidate genes for Mime (script 06)
load(file.path(OBJ_DIR, "bulk_degs_187.RData"))   # bulk_degs
write.csv(bulk_degs, file.path(result_path,"bulk_degs.csv"))
cand_new <- intersect(deg_sig$gene, bulk_degs$GeneSymbol)
cat("\nNEW candidates (malignant sc DEG ∩ bulk DEG):", length(cand_new), "\n")
print(sort(cand_new))


writeLines(cand_new, file.path(OBJ_DIR, "candidates_malignant_only_187tcga.txt"))
save(seu_mal, file = file.path(OBJ_DIR, "seu_mal_UPS_groups.RData"))

## 7. Malignant vs normal-like epithelial cells per sample
tab <- table(res$sample, res$malignant)
per_sample <- cbind(tab, Total = rowSums(tab),
                    Pct_malignant = round(100 * tab[, "TRUE"] / rowSums(tab), 1))
print(per_sample)
write.csv(per_sample, file.path(result_path, "malignant_cells_per_sample.csv"))

