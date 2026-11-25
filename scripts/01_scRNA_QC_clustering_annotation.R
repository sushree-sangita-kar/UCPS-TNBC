# 01  GSE176078 scRNA-seq: QC, clustering, Harmony integration, annotation,
#     epithelial sub-clustering
#
# Input : data/GSE176078/<sample>/count_matrix_sparse.mtx, _genes.tsv,
#         _barcodes.tsv, metadata.csv  (GEO GSE176078 supplementary files)
#         percent.mito is taken from the authors' metadata.csv
# Output: results/objects/seurat_all_annotated.RData   (object: seurat_harmony)
#         results/objects/epithelial_cells.RData       (object: epithelial_cells)
#         figures in results/01_scRNA_QC_clustering_annotation/

# Clear environment and setup
#rm(list = ls(all.names = TRUE)) # Clear all objects including hidden objects
#gc() # Free up memory and report memory usage

# Loading libraries
library(harmony)# load harmony before Seurat
library(Seurat)
library(scCustomize)
library(SingleCellExperiment)
library(sctransform)
library(patchwork)
library(pheatmap)

library(purrr)
library(readr)
library(stringr)
library(ragg)
library(tidyverse)
library(here)
library(msigdbr)
library(fgsea)
set.seed(123)

# Setting the paths
source(here::here("config.R"))
data_path  <- SC_DIR
plot_path  <- step_dir("01_scRNA_QC_clustering_annotation")
RData_path <- OBJ_DIR


# Load the data (one sub-folder per sample inside data/GSE176078)
dirs <- list.dirs(path = data_path, full.names = TRUE, recursive = FALSE)

# Load data
for (x in dirs) {
  dataset_name <- basename(x)
  mtx_file <- file.path(x, 'count_matrix_sparse.mtx')
  features_file <- file.path(x, 'count_matrix_genes.tsv')
  cells_file <- file.path(x, 'count_matrix_barcodes.tsv')
  metadata_file <- file.path(x, 'metadata.csv')
  
  if (file.exists(mtx_file) && file.exists(features_file) && file.exists(cells_file) && file.exists(metadata_file)) {
    count_matrix <- ReadMtx(
      mtx = mtx_file, features = features_file, feature.column = 1,cells = cells_file)
    metadata <- read.csv(metadata_file, row.names = 1)
    seurat_obj <- CreateSeuratObject(counts = count_matrix, meta.data = metadata)
    assign(dataset_name, seurat_obj)
  } else {
    cat("Skipping directory", dataset_name, "because one or more necessary files are missing.\n")
  }
}

#Merging only untreated TNBC samples(refer GEO176978 from GEO n supplimentary) into a object
# Merging the samples to create one seurat object
merged_seurat_objects <- merge( CID44991, y = c(CID44041, CID4465,CID4495, CID44971, CID4515),
                               project = 'brc')
unique(merged_seurat_objects@meta.data$orig.ident)
table(merged_seurat_objects$orig.ident)


# Filter poor quality cells
# Assess percent of mitochondrial counts in each cell 
# merged_seurat_objects$percent.mito <- PercentageFeatureSet(object = merged_seurat_objects,pattern = "^\(?i\)MT-"\)

# Violin plot to visualize QC matrices (before cell subsetting)
VlnPlot(object = merged_seurat_objects, pt.size = 0,
        features = c("nFeature_RNA", 
                     "nCount_RNA", 
                     "percent.mito"), 
        ncol = 3)

# Feature plots to visualize relationship between QC matrices
plot1 <- FeatureScatter(merged_seurat_objects,
                        feature1 = "nCount_RNA", 
                        feature2 = "percent.mito")
plot2 <- FeatureScatter(object = merged_seurat_objects,
                        feature1 = "nCount_RNA", 
                        feature2 = "nFeature_RNA")
plot1 + plot2

# Subset merged_seurat_objects_filtered (cell filtering)
merged_seurat_objects_filtered <- subset(merged_seurat_objects,
                                         subset = nFeature_RNA > 200 & nCount_RNA > 500 & percent.mito < 20)
merged_seurat_objects_filtered


# Violin plot to visualize QC matrices (after cell subsetting)
merged_seurat_objects_filtered <- JoinLayers(merged_seurat_objects_filtered)

md <- merged_seurat_objects_filtered@meta.data
# Rename orig.ident to Sample ID
md$`Sample ID` <- md$orig.ident


p1 <- ggplot(md, aes(x = `Sample ID`, y = nFeature_RNA, fill = `Sample ID`)) +
  geom_violin(trim = FALSE) + theme_classic() + 
  NoLegend() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
        axis.text.y = element_text(size = 16),
        axis.title = element_text(size = 17),
        plot.title = element_text(hjust = 0.5, size = 18))

p2 <- ggplot(md, aes(x = `Sample ID`, y = nCount_RNA, fill = `Sample ID`)) +
  geom_violin(trim = FALSE) + theme_classic() + 
  NoLegend() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
        axis.text.y = element_text(size = 16),
        axis.title = element_text(size = 17),
        plot.title = element_text(hjust = 0.5, size = 18))

p3 <- ggplot(md, aes(x = `Sample ID`, y = percent.mito, fill = `Sample ID`)) +
  geom_violin(trim = FALSE) + theme_classic() +
  NoLegend() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
        axis.text.y = element_text(size = 16),
        axis.title = element_text(size = 17),
        plot.title = element_text(hjust = 0.5, size = 18))

ggsave(filename = file.path(plot_path, "QC_violin_plots.png"), 
       plot = (p1 | p2 | p3), 
       width = 18, height = 5, dpi = 300)


# Visualizing the metadata
head(merged_seurat_objects_filtered@meta.data)



# Perform standard workflow steps 
merged_seurat_objects_filtered <- NormalizeData(object = merged_seurat_objects_filtered)
merged_seurat_objects_filtered <- FindVariableFeatures(object = merged_seurat_objects_filtered)
merged_seurat_objects_filtered <- ScaleData(object = merged_seurat_objects_filtered)
merged_seurat_objects_filtered <- RunPCA(object = merged_seurat_objects_filtered)

DimPlot_scCustom(merged_seurat_objects_filtered, reduction = "pca", label = TRUE, repel = TRUE)

# Choosing the best pc
graphics.off()
ElbowPlot(merged_seurat_objects_filtered, ndims=30, reduction="pca")
pct <- merged_seurat_objects_filtered$pca@stdev / sum(merged_seurat_objects_filtered$pca@stdev) * 100
cumu <- cumsum(pct)
co1 <- which(cumu > 90 & pct < 5)[1]
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
pcs <- min(co1, co2)
print(pcs)

# Cell clustering
merged_seurat_objects_filtered <- FindNeighbors(object = merged_seurat_objects_filtered, dims = 1:pcs)
merged_seurat_objects_filtered <- FindClusters(object = merged_seurat_objects_filtered,resolution = 0.2) 
merged_seurat_objects_filtered <- RunUMAP(object = merged_seurat_objects_filtered, dims = 1:pcs) 

# Check for batch effect (by coloring acc to sample_ID)
# UMAP visualization
png (file.path(plot_path,"all_cells_before_harmony.png"), width = 1200, height = 400)
DimPlot_scCustom(merged_seurat_objects_filtered, reduction = 'umap', label = TRUE, repel = TRUE, group.by = c( "celltype_major", "orig.ident"))
dev.off()



px<- DimPlot_scCustom(merged_seurat_objects_filtered, reduction = 'umap', group.by = c("seurat_clusters")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size = 23),
        plot.title = element_blank(),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 17))
ggsave(filename = file.path(plot_path, "seurat_clusters.png"), 
       plot =px, 
       width = 6, height = 4, dpi = 300)




p4<- DimPlot_scCustom(merged_seurat_objects_filtered, reduction = 'umap', group.by = c("orig.ident")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size = 23),
        plot.title = element_blank(),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 17))


ggsave(filename = file.path(plot_path, "Before batch effect correction.png"), 
       plot =p4, 
       width = 6, height = 4, dpi = 300)


# Batch correction using harmony
seurat_harmony <- harmony::RunHarmony(merged_seurat_objects_filtered, 
                             group.by.vars = "orig.ident")

# Re-run standard pipeline on Harmony-corrected data
seurat_harmony <- FindNeighbors(seurat_harmony, reduction = "harmony", dims = 1:20)
seurat_harmony <- FindClusters(seurat_harmony, resolution = 0.5)
seurat_harmony <- RunUMAP(seurat_harmony, reduction = "harmony", dims = 1:20)


png (file.path(plot_path,"all_cells_after_harmony.png"), width = 1200, height = 400)
DimPlot_scCustom(seurat_harmony, reduction = 'umap', label = TRUE, repel = TRUE, group.by = c( "celltype_major", "orig.ident"))
dev.off()



p6<- DimPlot_scCustom(seurat_harmony, reduction = 'umap', group.by = c("orig.ident"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size = 23),
        plot.title = element_blank(),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 17))

ggsave(filename = file.path(plot_path, "After batch effect correction.png"), 
       plot =p6, 
       width = 6, height = 4, dpi = 300)



py<- DimPlot_scCustom(seurat_harmony, reduction = 'umap', group.by = c("seurat_clusters"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size = 23),
        plot.title = element_blank(),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 17))

ggsave(filename = file.path(plot_path, "seurat_clusters_after_batch_effect_correction.png"), 
       plot =py, 
       width = 6, height = 5, dpi = 300)





# Find cluster-wise DEGs
# Check the unique clusters and number of cells per cluster
table(seurat_harmony@meta.data$seurat_clusters)


seurat_harmony[["RNA"]] <- JoinLayers(seurat_harmony[["RNA"]])

# Find DEGs for every cluster compared to all remaining cells, report only the positive ones
Idents(seurat_harmony) <- seurat_harmony@meta.data$seurat_clusters
markers <- FindAllMarkers(
  object = seurat_harmony,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25)

write.csv(markers, file.path(plot_path, "cluster_markers_all_cells.csv"), row.names = FALSE)

# Apply top_n within each cluster group to limit memory usage
top_markers <- markers %>%
  group_by(cluster) %>% 
  slice_max(avg_log2FC, n = 20)

# Manual annotation of cell types using known (canonical) cell markers
canonical_markers <- list(
  "Epithelial\nCells" = c("EPCAM", "CDH1"),
  "T Cells" = c("CD3E", "CD2", "CD3D"),
  "TAMs" = c("C1QC", "MRC1", "CCL18"),
  "CAFs" = c("PDGFRA", "FAP", "DCN"),
  "Pericytes" = c("MCAM", "PDGFRB", "RGS5"),
  "B Cells" = c("CD19", "CD79A", "MS4A1"),
  "Plasma\nCells" = c("IGHG1", "MZB1", "JCHAIN"),
  "Endothelial\nCells" = c("VWF", "CD34", "PECAM1"),
  "DCs" = c("IRF8", "CLECL1", "PTCRA")
)

Idents(seurat_harmony) <- "seurat_clusters"

p7<- DotPlot(seurat_harmony, features = canonical_markers) + RotatedAxis()
p7 <- p7 + ylab("Seurat clusters") + ggtitle("Cell type annotation using canonical markers")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size = 23),
        plot.title = element_blank(),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 17),
        strip.text = element_text(size = 16))

ggsave(file.path(plot_path, "Cell type annotation using canonical marker expression.svg"), 
       plot = p7, width = 20, height = 7, units = "in", dpi = 300)



head(Idents(seurat_harmony))

# Define the mapping between old cluster IDs and new cell type names
cluster_id_mapping <- c("1" = "Epithelial cells", 
                        "5" = "Epithelial cells", 
                        "4" = "Epithelial cells",
                        "13" = "Epithelial cells", 
                        "15" = "Epithelial cells",
                        "16" = "Epithelial cells",
                        "18" = "Epithelial cells",
                        "3" = "TAMs",
                        "17" = "DCs",
                        "11" = "Endothelial cells",
                        "7" = "B cells",
                        "10" = "T cells",
                        "9" = "T cells",
                        "2" = "T cells",
                        "0" = "T cells",
                        "6" = "CAFs",
                        "14" = "Pericytes",
                        "8" = "Plasma cells",
                        "12" = "Plasma cells")



# Rename the cluster IDs in  Seurat object
seurat_harmony <- RenameIdents(seurat_harmony, cluster_id_mapping)

# Add a new column named "annotation" to the metadata with the renamed cell types
seurat_harmony <- AddMetaData(seurat_harmony, 
                              metadata = as.data.frame(Idents(seurat_harmony)),
                              col.name = "annotation")

head(seurat_harmony@meta.data)
table(seurat_harmony$annotation, seurat_harmony$celltype_major)



# Visualize the clusters with the new names
png(file.path(plot_path,"manual_annotation.png"), width = 1200, height = 400)
DimPlot(seurat_harmony, reduction = "umap", group.by = "annotation" ,label = TRUE, pt.size = 0.5) 
dev.off()


# Fibroblast cluster is labelled "CAFs". The mapping above already uses "CAFs",
# so on a fresh run this block only checks the label.
ann <- seurat_harmony$annotation
old_levels <- if (is.factor(ann)) levels(ann) else sort(unique(ann))
new_levels <- sub("^Fibroblasts$", "CAFs", old_levels)          # same position = same colour
ann <- factor(sub("^Fibroblasts$", "CAFs", as.character(ann)), levels = new_levels)
seurat_harmony$annotation <- ann
Idents(seurat_harmony) <- "annotation"

table(seurat_harmony$annotation)

# All cell annotated
p8 <- DimPlot(seurat_harmony, reduction = "umap", group.by = "annotation", pt.size = 0.1) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size = 23),
        plot.title = element_blank(),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 17))

ggsave(file.path(plot_path, "All_cells_annotated.svg"), 
       plot = p8, width = 7, height = 4, units = "in", dpi = 300)



# Annotated object used by inferCNV (script 03) and CellChat (script 11)
save(seurat_harmony, file = file.path(RData_path, "seurat_all_annotated.RData"))


# Epithelial cells with harmony(batch effect removal)

# Step 1: Subset the Seurat object for epithelial cells
Idents(seurat_harmony) <- seurat_harmony$annotation
epithelial_cells <- subset(seurat_harmony, idents = c("Epithelial cells"))

# Step 2: Perform preprocessing
epithelial_cells <- NormalizeData(epithelial_cells)
epithelial_cells <- FindVariableFeatures(epithelial_cells, selection.method = "vst", nfeatures = 2000)
epithelial_cells <- ScaleData(epithelial_cells, features = rownames(epithelial_cells))

# Step 3: Run dimensionality reduction
epithelial_cells <- RunPCA(epithelial_cells, features = VariableFeatures(object = epithelial_cells))

# Choosing the best pc FIRST
ElbowPlot(epithelial_cells, ndims=30, reduction="pca")
pct <- epithelial_cells$pca@stdev / sum(epithelial_cells$pca@stdev) * 100
cumu <- cumsum(pct)
co1 <- which(cumu > 90 & pct < 5)[1]
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
pcs <- min(co1, co2)
print(pcs)

# Re-apply Harmony using calculated PCs
epithelial_cells <- harmony::RunHarmony(epithelial_cells, 
                                        group.by.vars = "orig.ident",
                                        dims.use = 1:pcs)  # Use calculated pcs

# Standard pipeline on harmony corrected data
epithelial_cells <- FindNeighbors(epithelial_cells, reduction = "harmony", dims = 1:pcs)
epithelial_cells <- FindClusters(epithelial_cells, resolution = 0.5)
epithelial_cells <- RunUMAP(epithelial_cells, reduction = "harmony", dims = 1:pcs)

# Step 5: Visualize the subclusters
png(file.path(plot_path,"Epithelial_cells_with_harmony.png"), width = 1200, height = 800)  # Fixed filename
DimPlot(epithelial_cells, reduction = "umap", group.by = c("annotation","orig.ident"), label = TRUE, pt.size = 0.1) 
dev.off()


# Epithelial object used for AUCell UPS scoring (script 02)
save(epithelial_cells, file = file.path(RData_path, "epithelial_cells.RData"))


# Load patchwork if not already loaded
library(patchwork)


# 1. Save Subclustering of epithelial cellsplot
p9a <- DimPlot(epithelial_cells, reduction = "umap", group.by = "annotation", pt.size = 0.1)  + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
                                                                                                                  axis.text.y = element_text(size = 20),
                                                                                                                  axis.title = element_text(size = 23),
                                                                                                                  plot.title = element_blank(),
                                                                                                                  legend.text = element_text(size = 16),
                                                                                                                  legend.title = element_text(size = 17))
p9b <- DimPlot(epithelial_cells, reduction = "umap", group.by = "seurat_clusters", pt.size = 0.1) + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
                                                                                                          axis.text.y = element_text(size = 20),
                                                                                                          axis.title = element_text(size = 23),
                                                                                                          plot.title = element_blank(),
                                                                                                          legend.text = element_text(size = 16),
                                                                                                          legend.title = element_text(size = 17))

p9 <- (p9a | p9b) 

ggsave(file.path(plot_path, "Subclustering of epithelial cells.png"), 
       plot = p9, width = 11, height = 4.5, units = "in", dpi = 300) 

# 2. Save Sample ID plot
p10 <- DimPlot(epithelial_cells, reduction = "umap", group.by = "orig.ident", pt.size = 0.1) + 
  ggtitle("Sample-specific clustering of epithelial cells") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
        
p10 <- p10 + theme(plot.title = element_text(hjust = 0.5),
                   plot.title.position = "plot")

ggsave(file.path(plot_path, "Epi_cells_after_batch_effect_correction.png"), 
       plot = p10, width = 6, height =4, units = "in", dpi = 300)



# THE END











