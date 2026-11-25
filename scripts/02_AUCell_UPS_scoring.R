# 02  UPS activity score per epithelial cell (AUCell)
#
# The scores used in this study are stored, one value per cell, in
#   resources/UPS_AUCell_scores_epithelial_cells.csv
# and script 05 reads them from there. This script shows how they were computed
#
# Input : results/objects/epithelial_cells.RData   (script 01)
#         resources/ups_genes.txt                  (UPS gene set, one gene per line)
# Output: results/02_AUCell_UPS/UPS_AUCell_computed.csv 

library(Seurat)
library(AUCell)
source(here::here("config.R"))
out_dir <- step_dir("02_AUCell_UPS")

load(file.path(OBJ_DIR, "epithelial_cells.RData"))   # epithelial_cells

## 1. UPS gene set present in the data
ups_all <- read.table(file.path(RES_DIR, "ups_genes.txt"))$V1
ups_genes_in_data <- intersect(ups_all, rownames(epithelial_cells))
cat("UPS genes in dataset:", length(ups_genes_in_data), "of", length(unique(ups_all)), "\n")
gene_sets <- list(UPS_signature = ups_genes_in_data)

## 2. AUCell on the normalised expression
expr_matrix <- as.matrix(GetAssayData(epithelial_cells, layer = "data", assay = "RNA"))

set.seed(42)   # makes a rerun of this script repeatable
cells_rankings <- AUCell_buildRankings(expr_matrix, plotStats = FALSE)
cells_AUC      <- AUCell_calcAUC(gene_sets, cells_rankings, nCores = 4)

aucell_scores <- as.numeric(getAUC(cells_AUC)["UPS_signature", ])
names(aucell_scores) <- colnames(expr_matrix)

recomputed <- data.frame(cell = names(aucell_scores),
                         UPS_AUCell = sprintf("%.17g", aucell_scores))
write.csv(recomputed, file.path(out_dir, "UPS_AUCell_computed.csv"),
          row.names = FALSE, quote = FALSE)



