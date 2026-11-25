# 04d  TCGA bulk DEGs (limma-voom): 187 tumours vs solid tissue normal
#      Tumour arm = the 187 patients used to train the Mime model.
#      Normal arm = TCGA solid-tissue-normal breast samples.
#
# Input : results/objects/TCGA_BRCA_data.RData                         (script 04a)
#         results/objects/187tcga_list_training_validation_Datasets.Rdata (script 04c)
# Output: results/objects/bulk_degs_187.RData   (object: bulk_degs; used by script 05)
#         results/04d_TCGA_bulk_DEG/  full table, volcano plot

# 0. Libraries
library(dplyr)
library(ggplot2)
library(edgeR)
library(limma)
library(SummarizedExperiment)
source(here::here("config.R"))
set.seed(1234)

RData_path  <- OBJ_DIR
plots_path  <- step_dir("04d_TCGA_bulk_DEG")

# 1. Load expression object + the MIME 187-sample cohort
load(file.path(OBJ_DIR, "TCGA_BRCA_data.RData"))   # -> `data` (SummarizedExperiment)

# The 187-patient list used in Mime:
load(file.path(OBJ_DIR, "187tcga_list_training_validation_Datasets.Rdata"))

# Patient IDs (12-char) that define the tumour arm.
mime_ids <- unique(substr(list_train_vali_Data$tcga$ID, 1, 12))
cat("MIME cohort unique patient IDs:", length(mime_ids), "\n")
stopifnot(length(mime_ids) == 187)   # hard guarantee we start from 187

# 2. PAM50 relabelling 
colData(data)$paper_BRCA_Subtype_PAM50[colData(data)$paper_BRCA_Subtype_PAM50 == "Normal" &
                                         colData(data)$tissue_type == "Tumor"] <- "Normal_like"
colData(data)$paper_BRCA_Subtype_PAM50[is.na(colData(data)$paper_BRCA_Subtype_PAM50) &
                                         colData(data)$tissue_type == "Normal"] <- "Normal"

# 3. Base filtering
data <- data[, !is.na(colData(data)$gender) &
               colData(data)$gender == "female" &
               !is.na(colData(data)$paper_BRCA_Subtype_PAM50) &
               (colData(data)$paper_BRCA_Subtype_PAM50 == "Normal" |
                  colData(data)$paper_BRCA_Subtype_PAM50 == "Basal")]

data <- data[, !is.na(colData(data)$vital_status) &
               ((colData(data)$vital_status == "Dead"  & !is.na(colData(data)$days_to_death)          & colData(data)$days_to_death > 0) |
                  (colData(data)$vital_status == "Alive" & !is.na(colData(data)$days_to_last_follow_up) & colData(data)$days_to_last_follow_up > 0))]

# 4. Restrict the tumour arm to the 187 Mime patients
# 12-char patient barcode for every remaining sample.
# Use colData$patient if present (already 12-char); else derive from barcode.
pat12 <- if (!is.null(colData(data)$patient)) {
  substr(as.character(colData(data)$patient), 1, 12)
} else {
  substr(colData(data)$barcode, 1, 12)
}

is_tumor  <- colData(data)$definition == "Primary solid Tumor"
is_normal <- colData(data)$definition == "Solid Tissue Normal"

# Tumour column kept iff its patient is one of the 187 MIME patients.
keep_tumor  <- is_tumor  & (pat12 %in% mime_ids)
keep_normal <- is_normal                                   # keep all normals as DE reference

data <- data[, keep_tumor | keep_normal]

# Collapse aliquot replicates: one tumour column per patient
pat12    <- if (!is.null(colData(data)$patient)) substr(as.character(colData(data)$patient),1,12) else substr(colData(data)$barcode,1,12)
is_tumor <- colData(data)$definition == "Primary solid Tumor"

# Flag a tumour as a duplicate ONLY if another TUMOUR of the same patient
# precedes it. Computing duplicated() within the tumour subset prevents a
# patient's matched NORMAL sample from making its own tumour look duplicated.
dup_tumor <- logical(ncol(data))
tumor_idx <- which(is_tumor)
dup_tumor[tumor_idx] <- duplicated(pat12[tumor_idx])   # within-tumour only
if (any(dup_tumor)) {
  cat("Dropping", sum(dup_tumor), "duplicate tumour aliquot(s)\n")
  data <- data[, !dup_tumor]
}

# Verify the tumour arm is exactly the 187 MIME patients
pat12     <- if (!is.null(colData(data)$patient)) substr(as.character(colData(data)$patient),1,12) else substr(colData(data)$barcode,1,12)
is_tumor  <- colData(data)$definition == "Primary solid Tumor"
tumor_pat <- pat12[is_tumor]

missing_from_data <- setdiff(mime_ids, tumor_pat)   # MIME patients with no expression column
if (length(missing_from_data) > 0) {
  warning("MIME patients absent from expression object: ",
          paste(missing_from_data, collapse = ", "))
}
cat("Tumour samples in DE:", sum(is_tumor),
    " | Normal samples in DE:", sum(!is_tumor), "\n")
stopifnot(length(unique(tumor_pat)) == length(tumor_pat))          # no dup patients
stopifnot(all(tumor_pat %in% mime_ids))                            # every tumour is a MIME patient
# If TCGA has an expression column for all 187, this must hold exactly:
stopifnot(sum(is_tumor) == length(setdiff(mime_ids, missing_from_data)))

table(colData(data)$definition, useNA = "always")

# 5. Gene ID -> Symbol mapping
# Gene symbols from the GENCODE annotation of the count object (rowData$gene_name),
# the same source used for log_tpm_matrix_median (script 04a). Same-symbol rows
# are collapsed by summing raw counts.
raw_count_mat <- assay(data, "unstranded")

gene_info <- rowData(data)
stopifnot(nrow(gene_info) == nrow(raw_count_mat))          # rowData aligns with counts
sym <- gene_info$gene_name

# Drop rows with no symbol or multi-gene (";") names, as for the log-TPM matrix.
ok  <- !is.na(sym) & sym != "" & !grepl(";", sym)
cat("Rows with usable gene_name:", sum(ok), "of", length(sym),
    sprintf(" (%.1f%% dropped)\n", 100*mean(!ok)))
raw_count_mat <- raw_count_mat[ok, ]
sym <- sym[ok]

# Collapse duplicate symbols by SUMMING raw counts across their Ensembl IDs.
n_dup_symbols <- sum(duplicated(sym))
cat("Symbols requiring sum-collapse (>=2 Ensembl IDs):",
    length(unique(sym[duplicated(sym)])),
    "| rows summed away:", n_dup_symbols, "\n")
raw_count_mat <- rowsum(raw_count_mat, group = sym)         # rownames become symbols

genes_df <- data.frame(
  GeneSymbol = rownames(raw_count_mat),
  row.names  = rownames(raw_count_mat)
)

stopifnot(!any(duplicated(rownames(raw_count_mat))))        # one row per symbol now

# 6. limma DEG (design built AFTER subsetting)
sample_info <- colData(data)
group  <- factor(sample_info$definition)
group  <- relevel(group, ref = "Solid Tissue Normal")
design <- model.matrix(~group)
stopifnot(nrow(design) == ncol(raw_count_mat))     # rows == samples

dge <- DGEList(counts = raw_count_mat, samples = colData(data), genes = genes_df)

keep <- filterByExpr(dge, design)                  # filter keyed to the FINAL design
dge  <- dge[keep, , keep.lib.sizes = FALSE]
rm(keep)

dge <- calcNormFactors(dge, method = "TMM")
v   <- voom(dge, design, plot = TRUE)

fit <- lmFit(v, design)
fit <- eBayes(fit)

topGenes <- topTable(fit, coef = 2, sort.by = "p", number = Inf)
# Select columns by name (topTable column positions depend on the annotation columns)
keep_cols <- intersect(c("GeneSymbol", "logFC", "AveExpr", "P.Value", "adj.P.Val", "B"),
                       colnames(topGenes))
topGenes  <- topGenes[, keep_cols]
stopifnot("adj.P.Val" %in% colnames(topGenes), "logFC" %in% colnames(topGenes))
save(topGenes, file = file.path(RData_path, "topGenes_187.RData"))
write.csv(topGenes, file.path(plots_path, "topGenes_187_full.csv"), row.names = FALSE)

# 7. Regulation calls + DEG export
topGenes$regulation <- "Not Significant"
topGenes$regulation[topGenes$logFC >  1 & topGenes$adj.P.Val < 0.05] <- "upregulated"
topGenes$regulation[topGenes$logFC < -1 & topGenes$adj.P.Val < 0.05] <- "downregulated"

bulk_degs <- topGenes[topGenes$regulation %in% c("upregulated", "downregulated"), ]
save(bulk_degs, file = file.path(RData_path, "bulk_degs_187.RData"))
write.csv(bulk_degs, file.path(plots_path, "bulk_degs_187.csv"), row.names = FALSE)

# 8. Volcano
p <- ggplot(topGenes, aes(x = logFC, y = -log10(adj.P.Val), color = regulation)) +
  geom_point(size = 0.8) +
  scale_color_manual(values = c("upregulated"     = "blue",
                                "downregulated"   = "red",
                                "Not Significant" = "gray")) +
  labs(title = "TCGA DEGs Volcano Plot (187 MIME cohort)",
       x = "log2 Fold Change", y = "-log10 Adjusted P-value") +
  theme_minimal() +
  theme(legend.position = "right",
        plot.title  = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.title  = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 10, face = "bold"))
ggsave(file.path(plots_path, "Volcano_Plot_of_TCGA_DEGs_187.png"),
       plot = p, width = 8, height = 6, dpi = 120)
