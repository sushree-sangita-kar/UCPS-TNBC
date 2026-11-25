## 03  inferCNV per sample: epithelial cells (observations) vs all non-epithelial
##     cell types (reference), then malignant calls
##
## The cell labels used in the published inferCNV run are stored in
## resources/inferCNV_cell_annotation.csv. The malignant calls
## used in the study are stored in resources/malignant_classification_anchored.csv;
## script 05 reads the stored calls.
##
## Input : results/objects/seurat_all_annotated.RData   (script 01)
##         resources/inferCNV_cell_annotation.csv
##         data/reference/hg38_gencode_v27.txt          (inferCNV gene order file)
## Output: results/03_inferCNV/<sample>/output/          (inferCNV run files)
##         results/objects/infercnv/infercnv_<sample>.RData
##         results/03_inferCNV/malignant_classification_recomputed.csv

library(Seurat)
library(infercnv)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)

# Set options
options(scipen = 100)
set.seed(1234)

## 1. CONFIGURATION

# Paths
source(here::here("config.R"))
sc_data         <- file.path(OBJ_DIR, "seurat_all_annotated.RData")
gene_order_file <- file.path(REF_DIR, "hg38_gencode_v27.txt")
base_dir        <- step_dir("03_inferCNV")
rdata_dir       <- file.path(OBJ_DIR, "infercnv")
dir.create(rdata_dir, recursive = TRUE, showWarnings = FALSE)

# Cell labels of the published inferCNV run
apply_infercnv_labels <- function(obj) {
  lab <- read.csv(file.path(RES_DIR, "inferCNV_cell_annotation.csv"))
  stopifnot(setequal(lab$cell, colnames(obj)))
  obj$annotation <- lab$annotation[match(colnames(obj), lab$cell)]
  obj
}

# Parameters
epi_label <- "Epithelial cells"
samples <- c("CID44041", "CID4465", "CID4495", "CID44971", "CID44991", "CID4515")

## 2. LOAD DATA

cat("Loading Seurat object...\n")
load(sc_data)
seurat_harmony <- apply_infercnv_labels(seurat_harmony)

cat("Total cells:", ncol(seurat_harmony), "\n")
cat("Samples:", paste(samples, collapse = ", "), "\n\n")

## 3. RUN INFERCNV FOR EACH SAMPLE

for (sample_name in samples) {
  
  cat("\n", rep("=", 60), "\n", sep = "")
  cat("Processing Sample:", sample_name, "\n")
  cat(rep("=", 60), "\n\n", sep = "")
  
  # Create sample-specific directories
  sample_dir <- file.path(base_dir, sample_name)
  dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Subset cells for this sample
  seurat_sample <- subset(seurat_harmony, subset = orig.ident == sample_name)
  n_cells <- ncol(seurat_sample)
  
  cat("Total cells in", sample_name, ":", n_cells, "\n")
  
  # Get cell annotations
  cell_annotations <- seurat_sample@meta.data$annotation
  names(cell_annotations) <- colnames(seurat_sample)
  
  # Remove NA annotations
  valid_cells <- names(cell_annotations)[!is.na(cell_annotations)]
  cell_annotations <- cell_annotations[valid_cells]
  
  # Check cell type distribution
  cell_type_table <- table(cell_annotations)
  cat("\nCell type distribution:\n")
  print(cell_type_table)
  
  # Define reference groups (non-epithelial)
  ref_groups <- setdiff(unique(cell_annotations), epi_label)
  n_epithelial <- sum(cell_annotations == epi_label)
  n_reference <- sum(cell_annotations != epi_label)
  
  cat("\nEpithelial cells:", n_epithelial)
  cat("\nReference cells:", n_reference)
  cat("\nReference groups:", length(ref_groups), "\n")
  
  # Skip if insufficient cells
  if (length(ref_groups) == 0 || n_epithelial < 20) {
    cat("\nWARNING: Insufficient cells for", sample_name, "- skipping\n")
    next
  }
  
  tryCatch({
    
    ## 3.1 Prepare expression matrix
    
    # Extract raw counts
    counts <- GetAssayData(seurat_sample, assay = "RNA", slot = "counts")
    
    # Convert sparse to dense if needed
    if (inherits(counts, "dgCMatrix")) {
      counts <- as.matrix(counts)
    }
    
    # Filter to valid cells only
    counts <- counts[, valid_cells, drop = FALSE]
    
    cat("\nExpression matrix:", nrow(counts), "genes x", ncol(counts), "cells\n")
    
    ## 3.2 Create annotation file
    
    ann_df <- data.frame(
      cell = valid_cells,
      cell_type = cell_annotations[valid_cells],
      stringsAsFactors = FALSE
    )
    
    
    # Remove cells from cell types having less than 2 cells
    cell_type_counts <- table(ann_df$cell_type)
    valid_ref_types <- names(cell_type_counts)[cell_type_counts >= 2]
    ref_groups <- setdiff(valid_ref_types, epi_label)
    
    ann_df <- ann_df[ann_df$cell_type %in% c(epi_label, ref_groups), ]
    counts <- counts[, ann_df$cell, drop = FALSE]
    
    # Write annotation file (NO HEADER for inferCNV)
    ann_file <- file.path(sample_dir, "annotations.txt")
    write.table(ann_df, ann_file, sep = "\t", quote = FALSE, 
                row.names = FALSE, col.names = FALSE)
    
    ## 3.3 Create inferCNV object
    
    cat("\nCreating inferCNV object...\n")
    
    infercnv_obj <- CreateInfercnvObject(
      raw_counts_matrix = counts,
      annotations_file = ann_file,
      delim = "\t",
      gene_order_file = gene_order_file,
      ref_group_names = ref_groups
    )
    
    ## 3.4 Run inferCNV
    
    out_dir <- file.path(sample_dir, "output")
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    
    cat("Running inferCNV analysis...\n")
    cat("This may take 30-60 minutes...\n")
    
    infercnv_obj_run <- infercnv::run(
      infercnv_obj,
      cutoff = 0.1,
      out_dir = out_dir,
      cluster_by_groups = TRUE,
      denoise = TRUE,
      HMM = TRUE,
      HMM_type = "i6",
      analysis_mode = "samples",
      tumor_subcluster_partition_method = "leiden",
      num_threads = 1,
      no_plot = FALSE,
      output_format = "pdf",
      save_rds = TRUE,
      save_final_rds = TRUE,
      resume_mode = TRUE,
      png_res = 300
    )
    
    ## 3.5 Save the inferCNV object (malignant calls are made in section 4)
    if (file.exists(file.path(out_dir, "run.final.infercnv_obj"))) {
      save(infercnv_obj_run,
           file = file.path(rdata_dir, paste0("infercnv_", sample_name, ".RData")))
    }
    
    cat("\n✅ Completed", sample_name, "\n")
    
  }, error = function(e) {
    cat("\n❌ ERROR in", sample_name, ":", conditionMessage(e), "\n")
  })
  
  # Clean memory
  gc()
}


## 4. Malignant call per sample (reference-anchored threshold)
all_res <- list()
for (s in samples) {
  load(file.path(rdata_dir, paste0("infercnv_", s, ".RData")))   # -> infercnv_obj_run
  obj  <- infercnv_obj_run
  expr <- obj@expr.data                       # genes x ALL cells (denoised)
  
  ## reference (non-epithelial) vs observation (epithelial) cells, from the object
  ref_cells <- colnames(expr)[ unlist(obj@reference_grouped_cell_indices) ]
  epi_cells <- colnames(expr)[ unlist(obj@observation_grouped_cell_indices) ]
  
  ## CNV signal = mean squared deviation from the diploid baseline (1)
  sig <- colMeans((expr - 1)^2)
  
  ## Reference-anchored threshold: 95th percentile of the reference cells
  thr <- quantile(sig[ref_cells], 0.95)
  
  df <- data.frame(cell = epi_cells, sample = s,
                   cnv_signal = sig[epi_cells],
                   malignant = sig[epi_cells] > thr)
  all_res[[s]] <- df
  cat(sprintf("%s: malignant %d / %d epithelial | ref-95%% thr=%.4f | median epi sig=%.4f\n",
              s, sum(df$malignant), nrow(df), thr, median(sig[epi_cells])))
}
res <- do.call(rbind, all_res)
cat("\nOVERALL:", sum(res$malignant), "/", nrow(res), "epithelial cells malignant\n")
write.csv(res, file.path(base_dir, "malignant_classification_recomputed.csv"), row.names = FALSE)
