# 09  TCGA (n = 187) immune microenvironment by UCPS risk group
#     ESTIMATE -> CIBERSORTx (LM22) -> immune-checkpoint gene expression
#
# Input : results/objects/187tcga_list_training_validation_Datasets.Rdata (script 04c)
#         results/objects/tcga_risk_malignant_only_tcga187.RData         (script 06)
#         results/09_immune_ESTIMATE_CIBERSORT/CIBERSORTx_Job6_Results.csv
#           CIBERSORTx is run on the web portal (cibersortx.stanford.edu) with the
#           mixture file this script writes (TCGA_for_CIBERSORT.txt); downloaded the
#           results CSV into the folder above before running section 3.
# Output: results/09_immune_ESTIMATE_CIBERSORT/

# Packages
library(limma)
library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(survival)
library(survminer)
library(GSVA)
library(pheatmap)
library(corrplot)
library(dplyr)
library(tidyr)
library(readr)
library(tibble)
library(estimate)

# Paths
source(here::here("config.R"))
load(file.path(OBJ_DIR, "187tcga_list_training_validation_Datasets.Rdata"))
load(file.path(OBJ_DIR, "tcga_risk_malignant_only_tcga187.RData"))   # tcga_risk
RData_path <- OBJ_DIR
plot_path  <- step_dir("09_immune_ESTIMATE_CIBERSORT")
file_path  <- plot_path

# 1. Prepare data
tcga_expr <- list_train_vali_Data[["tcga"]]

if (any(nchar(tcga_expr$ID) > 12)) {
  sample_type <- substr(tcga_expr$ID, 14, 15)
  cat("Sample-type codes present:\n"); print(table(sample_type))
  tcga_expr <- tcga_expr[sample_type == "01", ]          # primary solid tumor
}

# Trim aliquot barcode -> patient barcode
tcga_expr$ID <- substr(tcga_expr$ID, 1, 12)

# Duplicates should now reflect genuine technical replicates only
cat("Duplicates after trim:", sum(duplicated(tcga_expr$ID)), "\n")
#tcga_expr <- tcga_expr[!duplicated(tcga_expr$ID), ]

# Verify overlap before merging
cat("Overlap with tcga_risk:", length(intersect(tcga_expr$ID, tcga_risk$ID)), "\n")

# Merge risk group / risk score
tcga_expr <- merge(tcga_expr,
                   tcga_risk[, c("ID", "risk_group", "RS")],
                   by = "ID")

cat("Final merged rows:", nrow(tcga_expr), "\n")
print(table(tcga_expr$risk_group))

# Build gene x sample expression matrix (log2(TPM+1) space)
meta_cols <- c("ID", "OS.time", "OS", "risk_group", "RS")
expr_df   <- tcga_expr[, !colnames(tcga_expr) %in% meta_cols]

# Guard against a stray non-numeric column silently coercing the
# whole matrix to character on t().
stopifnot(all(vapply(expr_df, is.numeric, logical(1))))

expr_mat <- as.matrix(expr_df)
rownames(expr_mat) <- tcga_expr$ID
expr_mat <- t(expr_mat)                      # genes x samples, log2(TPM+1)

save(expr_mat, file = file.path(RData_path, "expr_mat_tcga.RData"))   # log2(TPM+1), genes x samples

# 2. IMMUNE INFILTRATION — ESTIMATE 
estimate_input <- expr_mat
estimate_file  <- file.path(file_path, "estimate_input.txt")
write.table(estimate_input, estimate_file, sep = "\t", quote = FALSE,
            col.names = NA)               # explicit corner cell -> clean header

filterCommonGenes(input.f = estimate_file,
                  output.f = file.path(file_path, "estimate_filtered.gct"),
                  id = "GeneSymbol")

estimateScore(input.ds  = file.path(file_path, "estimate_filtered.gct"),
              output.ds = file.path(file_path, "estimate_scores.gct"),
              platform  = "illumina")       # RNA-seq; no affymetrix purity model

# Read ESTIMATE scores
estimate_scores <- read.table(file.path(file_path, "estimate_scores.gct"),
                              skip = 2, header = TRUE, row.names = 1)
estimate_scores <- t(estimate_scores[, -1])          # drop Description; samples x scores
estimate_scores <- as.data.frame(estimate_scores)
rownames(estimate_scores) <- gsub("\\.", "-", rownames(estimate_scores))
estimate_scores$ID <- substr(rownames(estimate_scores), 1, 12)   # patient barcode

# Merge with risk groups
estimate_df <- merge(estimate_scores, tcga_risk[, c("ID", "risk_group")], by = "ID")

# Compare ESTIMATE scores between groups
estimate_melt <- reshape2::melt(estimate_df, id.vars = c("ID", "risk_group"),
                                variable.name = "Score_Type", value.name = "Score")

p_est <- ggboxplot(estimate_melt, x = "risk_group", y = "Score",
                   fill = "risk_group", palette = c("#C0392B", "#1F77B4"),
                   facet.by = "Score_Type", scales = "free_y") +
  stat_compare_means(method = "wilcox.test", label = "p.format", size = 7) +
  labs(x = "Risk Group", y = "Score") +
  theme(legend.position = "bottom",
        axis.title = element_text(size = 18, face = "bold"),
        axis.text  = element_text(size = 16),
        strip.text = element_text(size = 16, face = "bold"),
        legend.text  = element_text(size = 15),
        legend.title = element_text(size = 16, face = "bold"))

ggsave(file.path(plot_path, "estimate_scores_comparison.png"),
       p_est, width = 10, height = 8, dpi = 300)

# 3. IMMUNE INFILTRATION — CIBERSORT (LM22)
# CIBERSORT requires LINEAR (non-log) mixture input. expr_mat is
# log2(TPM+1); back-transform to TPM before writing the CIBERSORT mixture file.
cibersort_lin  <- 2^expr_mat - 1
cibersort_input <- as.data.frame(cibersort_lin)
cibersort_input$GeneSymbol <- rownames(cibersort_input)
cibersort_input <- cibersort_input[, c("GeneSymbol",
                                       setdiff(colnames(cibersort_input), "GeneSymbol"))]

# Write to file_path, consistent with the rest of the pipeline
write.table(cibersort_input,
            file = file.path(file_path, "TCGA_for_CIBERSORT.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Read CIBERSORTx result (run externally on the linear-space mixture)
res <- read_csv(file.path(file_path, "CIBERSORTx_Job6_Results.csv"))

# Deconvolution-confidence filter; report drop counts per group
n_before <- nrow(res)
res <- res[res$`P-value` < 0.05, ]
cat("CIBERSORT samples passing p<0.05:", nrow(res), "/", n_before, "\n")

# Keep Mixture + all 22 LM22 cell types by NAME (drop the 3 stat
# columns explicitly) so column position can never drop a cell type.
immune_fraction <- res[, !colnames(res) %in% c("P-value", "Correlation", "RMSE")]

# Rename Mixture -> ID and reduce to patient barcode
colnames(immune_fraction)[colnames(immune_fraction) == "Mixture"] <- "ID"
immune_fraction$ID <- gsub("\\.", "-", immune_fraction$ID)
immune_fraction$ID <- substr(immune_fraction$ID, 1, 12)

# Merge with risk group
immune_merged <- merge(immune_fraction,
                       tcga_risk[, c("ID", "risk_group")],
                       by = "ID")

# Report per-group n so uneven filtering is visible
cat("CIBERSORT n per risk group:\n"); print(table(immune_merged$risk_group))

# Per-cell-type stats with BH correction
immune_long <- immune_merged %>%
  pivot_longer(cols = -c(ID, risk_group),
               names_to = "CellType", values_to = "Fraction")

cibersort_stats <- immune_long %>%
  group_by(CellType) %>%
  summarise(p_value = wilcox.test(Fraction ~ risk_group)$p.value,
            .groups = "drop") %>%
  mutate(adj_p = p.adjust(p_value, method = "BH")) %>%
  arrange(adj_p)

print(cibersort_stats)
write.csv(cibersort_stats,
          file.path(file_path, "cibersort_wilcox_BH.csv"), row.names = FALSE)

# Boxplot annotated with BH-ADJUSTED p
# Build significance labels from adjusted p and place them with stat_pvalue_manual
sig_lab <- cibersort_stats %>%
  mutate(group1 = "High", group2 = "Low",
         p.signif = cut(adj_p,
                        breaks = c(-Inf, 1e-4, 1e-3, 1e-2, 5e-2, Inf),
                        labels = c("****", "***", "**", "*", "ns"))) %>%
  filter(p.signif != "ns")

# y position for each label = just above the max fraction of that cell type
ypos <- immune_long %>% group_by(CellType) %>%
  summarise(y.position = max(Fraction) * 1.05, .groups = "drop")
sig_lab <- left_join(sig_lab, ypos, by = "CellType")

p_cib <- ggplot(immune_long, aes(x = CellType, y = Fraction, fill = risk_group)) +
  geom_boxplot(width = 0.7, position = position_dodge(0.8), outlier.shape = NA) +
  scale_fill_manual(values = c("High" = "#E06B5A", "Low" = "#5B9BD5")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.12))) +
  theme_classic(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 13),
        axis.text.y = element_text(size = 13),
        axis.title.y = element_text(size = 16),
        legend.position = "top",
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 15)) +
  labs(x = NULL, y = "Cell proportion", fill = "Risk group")

if (nrow(sig_lab) > 0) {
  p_cib <- p_cib +
    stat_pvalue_manual(sig_lab, x = "CellType", label = "p.signif",
                       y.position = "y.position", size = 6)
}

# Save once, one set of dimensions.
ggsave(file.path(plot_path, "cibersort_risk_group.png"),
       p_cib, width = 9, height = 5, dpi = 300)

# 4. Immune-checkpoint gene expression (log2(TPM+1) space is correct here)
checkpoint_genes <- c(
  "PDCD1", "CD274", "PDCD1LG2", "CTLA4", "LAG3", "HAVCR2", "TIGIT", "BTLA",
  "VSIR", "IDO1", "CD276", "TNFRSF9", "ICOS", "TNFRSF18", "CD27", "CD28"
)

available_genes <- intersect(checkpoint_genes, rownames(expr_mat))
cat("Found:", length(available_genes), "/", length(checkpoint_genes), "\n")
cat("Missing:", setdiff(checkpoint_genes, available_genes), "\n")

checkpoint_expr <- expr_mat[available_genes, , drop = FALSE] %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Gene") %>%
  pivot_longer(cols = -Gene, names_to = "ID", values_to = "Expression")

checkpoint_expr$ID <- gsub("\\.", "-", checkpoint_expr$ID)
checkpoint_expr$ID <- substr(checkpoint_expr$ID, 1, 12)

checkpoint_expr <- merge(checkpoint_expr,
                         tcga_risk[, c("ID", "risk_group")], by = "ID")
checkpoint_expr$risk_group <- factor(checkpoint_expr$risk_group,
                                     levels = c("High", "Low"))

# Multiple-testing correction across checkpoint genes (was absent)
checkpoint_stats <- checkpoint_expr %>%
  group_by(Gene) %>%
  summarise(p_value = wilcox.test(Expression ~ risk_group)$p.value,
            .groups = "drop") %>%
  mutate(adj_p = p.adjust(p_value, method = "BH")) %>%
  arrange(adj_p)
print(checkpoint_stats)
write.csv(checkpoint_stats,
          file.path(file_path, "checkpoint_wilcox_BH.csv"), row.names = FALSE)

ck_lab <- checkpoint_stats %>%
  mutate(group1 = "High", group2 = "Low",
         p.signif = cut(adj_p,
                        breaks = c(-Inf, 1e-4, 1e-3, 1e-2, 5e-2, Inf),
                        labels = c("****", "***", "**", "*", "ns"))) %>%
  filter(p.signif != "ns")
ck_ypos <- checkpoint_expr %>% group_by(Gene) %>%
  summarise(y.position = max(Expression) * 1.05, .groups = "drop")
ck_lab <- left_join(ck_lab, ck_ypos, by = "Gene")

p_checkpoint <- ggplot(checkpoint_expr,
                       aes(x = Gene, y = Expression, fill = risk_group)) +
  geom_boxplot(width = 0.7, position = position_dodge(0.8), outlier.shape = NA) +
  scale_fill_manual(values = c("High" = "#E06B5A", "Low" = "#5B9BD5")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.12))) +
  labs(x = NULL, y = "Expression (log2(TPM + 1))", fill = "Risk group") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 13, face = "italic"),
        axis.text.y = element_text(size = 13),
        axis.title.y = element_text(size = 15, face = "bold"),
        legend.position = "top",
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 14, face = "bold"))

if (nrow(ck_lab) > 0) {
  p_checkpoint <- p_checkpoint +
    stat_pvalue_manual(ck_lab, x = "Gene", label = "p.signif",
                       y.position = "y.position", size = 6)
}

ggsave(file.path(plot_path, "checkpoint_expression_boxplot.png"),
       p_checkpoint, width = 12, height = 6, dpi = 300)

