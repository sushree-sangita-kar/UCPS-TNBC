# 10  Drug-response prediction (oncoPredict, GDSC2) by UCPS risk group, TCGA n = 187
#     GDSC2 -> calcPhenotype -> predicted ln(IC50), High vs Low risk (Wilcoxon)
#
# Input : data/oncoPredict/Training Data/  GDSC2 expression + response files
#           (download from https://osf.io/c6tfx/)
#         results/objects/187tcga_list_training_validation_Datasets.Rdata (script 04c)
#         results/objects/tcga_risk_malignant_only_tcga187.RData         (script 06)
# Output: results/10_oncoPredict/

suppressPackageStartupMessages({
  library(oncoPredict)
  library(ggpubr)
  library(reshape2)
  library(dplyr)
  library(ggplot2)
})

## CONFIG
source(here::here("config.R"))
data_dir   <- file.path(DATA_DIR, "oncoPredict", "Training Data")
expr_rdata <- file.path(OBJ_DIR, "187tcga_list_training_validation_Datasets.Rdata")
risk_file  <- file.path(OBJ_DIR, "tcga_risk_malignant_only_tcga187.RData")
out_dir    <- step_dir("10_oncoPredict")
setwd(out_dir)   # calcPhenotype writes calcPhenotype_Output/ into the working directory

## 1. GDSC2 training data (download: https://osf.io/c6tfx/)
trainingExprData <- readRDS(file.path(data_dir, "GDSC2_Expr (RMA Normalized and Log Transformed).rds"))
trainingPtype    <- readRDS(file.path(data_dir, "GDSC2_Res.rds"))


trainingPtype <- exp(trainingPtype)

## 2. Test expression (TCGA) — full matrix
load(expr_rdata)                                     # list_train_vali_Data
load(risk_file)                                      # tcga_risk (ID, risk_group, RS)
stopifnot(all(c("ID","risk_group","RS") %in% colnames(tcga_risk)))

tcga_expr    <- list_train_vali_Data[["tcga"]]
tcga_expr$ID <- substr(tcga_expr$ID, 1, 12)          # patient-level ID (primary tumor)
tcga_expr    <- tcga_expr[!duplicated(tcga_expr$ID), ]
gene_cols    <- setdiff(colnames(tcga_expr), c("ID","OS.time","OS","risk_group","RS"))
expr_mat     <- t(as.matrix(tcga_expr[, gene_cols]))
storage.mode(expr_mat) <- "double"
colnames(expr_mat) <- tcga_expr$ID

# Clean duplicate artifacts from symbols
base_sym        <- sub("\\.[0-9]+$", "", gene_cols)
is_dup_artifact <- (base_sym != gene_cols) & (base_sym %in% gene_cols)
rownames(expr_mat) <- ifelse(is_dup_artifact, base_sym, gene_cols)
rv       <- rowSums((expr_mat - rowMeans(expr_mat, na.rm = TRUE))^2, na.rm = TRUE)
ord      <- order(rv, decreasing = TRUE)
expr_mat <- expr_mat[ord, , drop = FALSE][!duplicated(rownames(expr_mat)[ord]), , drop = FALSE]

testExprData <- expr_mat
if (max(testExprData, na.rm = TRUE) > 50) stop("test matrix looks LINEAR-scale: expect log2(x+1)")

cat("Training expr:", dim(trainingExprData), "| Training drug:", dim(trainingPtype), "\n")
cat("Test expr    :", dim(testExprData), "(full gene set, 187 patients)\n")
n_overlap <- length(intersect(rownames(trainingExprData), rownames(testExprData)))
cat(sprintf("Gene overlap (train vs test): %d\n", n_overlap))
stopifnot(n_overlap > 5000)

## 3. Predict drug response
drug_pred <- calcPhenotype(
  trainingExprData = trainingExprData,
  trainingPtype    = trainingPtype,
  testExprData     = testExprData,
  batchCorrect     = "eb",
  powerTransformPhenotype = TRUE,
  removeLowVaryingGenes    = 0.2,
  removeLowVaringGenesFrom = "homogenizeData",
  minNumSamples    = 10,
  printOutput      = TRUE,
  folder           = TRUE,          
  selection        = 1
)

drug_pred <- read.csv("calcPhenotype_Output/DrugPredictions.csv",
                      row.names = 1, check.names = FALSE)
cat("Predicted:", nrow(drug_pred), "samples x", ncol(drug_pred), "drugs\n")

# Predicted IC50 -> ln(IC50); lower values = higher predicted sensitivity
drug_pred <- log(drug_pred)

## 4. Merge with risk groups (Global Discovery Screen)
# Retain all predicted columns instead of filtering for a restricted panel
ids <- rownames(drug_pred)
if (all(nchar(ids) >= 15)) {
  keep <- substr(ids, 14, 15) == "01"
  cat("Full barcodes detected; primary-tumor kept:", sum(keep), "of", length(keep), "\n")
  drug_pred <- drug_pred[keep, , drop = FALSE]
  drug_pred$ID <- substr(rownames(drug_pred), 1, 12)
  drug_pred <- drug_pred[!duplicated(drug_pred$ID), ]
} else {
  drug_pred$ID <- rownames(drug_pred)
}

drug_df <- merge(drug_pred, tcga_risk[, c("ID","risk_group","RS")], by = "ID")
drug_df$risk_group <- factor(drug_df$risk_group, levels = c("Low","High"))
cat("Merged rows:", nrow(drug_df), "\n")
print(table(drug_df$risk_group))

## 5. Stats: Wilcoxon (High vs Low), BH across all 198 drug models

drug_melt <- reshape2::melt(drug_df, id.vars = c("ID","risk_group","RS"),
                            variable.name = "Drug_ID", value.name = "IC50")
drug_melt$Drug_ID <- as.character(drug_melt$Drug_ID)

lab_map <- data.frame(Drug_ID = unique(drug_melt$Drug_ID))
lab_map$base <- sub("_\\d+$", "", lab_map$Drug_ID)
is_dup <- lab_map$base %in% lab_map$base[duplicated(lab_map$base)]
lab_map$Drug <- ifelse(is_dup,
                       paste0(lab_map$base, " (", sub("^.*_", "", lab_map$Drug_ID), ")"),
                       lab_map$base)
cat("Compounds screened twice:", sum(duplicated(lab_map$base)), "\n")
stopifnot(!any(duplicated(lab_map$Drug)))
drug_melt <- left_join(drug_melt, lab_map[, c("Drug_ID","Drug")], by = "Drug_ID")

drug_stats <- drug_melt %>%
  group_by(Drug_ID, Drug) %>%
  summarise(
    p_value     = wilcox.test(IC50 ~ risk_group)$p.value,
    median_high = median(IC50[risk_group == "High"], na.rm = TRUE),
    median_low  = median(IC50[risk_group == "Low"], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    delta_median   = median_high - median_low,
    p_adj          = p.adjust(p_value, method = "BH"),
    more_sensitive = ifelse(delta_median < 0, "High-risk", "Low-risk")
  ) %>%
  arrange(p_value)

# Numbers reported in the Results: drugs tested, P < 0.05, BH q < 0.05, and the
# direction among the BH-significant drugs (more_sensitive = group with the
# lower median predicted ln(IC50))
cat("Drugs tested:", nrow(drug_stats),
    "| P < 0.05:", sum(drug_stats$p_value < 0.05),
    "| BH q < 0.05:", sum(drug_stats$p_adj < 0.05), "\n")
print(table(more_sensitive = drug_stats$more_sensitive[drug_stats$p_adj < 0.05]))

write.csv(drug_stats, file.path(out_dir, "drug_stats_all_drugs.csv"), row.names = FALSE)
write.csv(drug_stats[drug_stats$p_adj < 0.05, ],
          file.path(out_dir, "drug_stats_BH_significant.csv"), row.names = FALSE)

## 6. Top 10 drugs (lowest Wilcoxon P, all BH q < 0.05): boxplots
candidate_drugs <- drug_stats %>% filter(p_adj < 0.05) %>% arrange(p_value)
top10 <- head(candidate_drugs, 10)
top_paper_drugs <- top10$Drug
print(as.data.frame(top10))

plot_data <- drug_melt %>% filter(Drug %in% top_paper_drugs)

# Panels in order of significance
plot_data$Drug <- factor(plot_data$Drug, levels = top_paper_drugs)

# Formatting function for nominal p-values
fmt_p <- function(p) ifelse(p < 1e-3, formatC(p, format = "e", digits = 2),
                            formatC(p, format = "f", digits = 3))

# Build annotation table using p_value to match  plot labels
ann <- plot_data %>% 
  group_by(Drug) %>%
  summarise(y.position = max(IC50, na.rm = TRUE) * 1.05, .groups = "drop") %>%
  left_join(drug_stats %>% transmute(Drug, p_value), by = "Drug") %>%
  mutate(group1 = "Low", group2 = "High", label = paste0("p = ", fmt_p(p_value)))

# Generate the boxplot grid 
p_box <- ggboxplot(plot_data, x = "risk_group", y = "IC50",
                   fill = "risk_group", palette = c(Low = "#2166AC", High = "#B2182B"),
                   facet.by = "Drug", scales = "free_y", ncol = 5) +
  stat_pvalue_manual(ann, label = "label", tip.length = 0.015, size = 3.5) +
  labs(y = "Predicted ln(IC50)", x = "Risk Group", fill = "Risk Group") +
  theme_bw() +
  theme(
    legend.position = "none", 
    strip.background = element_rect(fill="gray95"),
    strip.text = element_text(face = "bold", size = 14),
    axis.text.x = element_text(size = 13)
  )

ggsave(filename = file.path(out_dir, "top10_drugs_IC50_by_risk_group.png"), 
       plot = p_box, width = 10, height = 7)


cat("Analysis complete. Check files in:", out_dir, "\n")


## 7. Spearman correlation of predicted ln(IC50) with the risk score
##         (same 10 drugs)
# 1. Predicted ln(IC50) of the same 10 drug models (by GDSC2 drug ID)
verified_clean_names <- top10$Drug
drug_df_clean <- drug_df[, c("ID", "RS", "risk_group", top10$Drug_ID), drop = FALSE]
colnames(drug_df_clean)[-(1:3)] <- verified_clean_names

# 2. Extract the raw gene expression data for your 6 risk genes safely
risk_genes      <- c("RARRES1", "KRT16", "ASPHD1", "EXOSC4", "HELLS", "PUF60")
available_genes <- intersect(risk_genes, rownames(testExprData))
gene_expr_subset <- as.data.frame(t(testExprData[available_genes, , drop = FALSE]))
gene_expr_subset$ID <- rownames(gene_expr_subset)

# 3. Merge gene expression, risk scores, and predicted drug data into a clean matrix
matrix_df <- merge(drug_df_clean, gene_expr_subset, by = "ID")

scatter_data <- reshape2::melt(matrix_df, id.vars = c("ID", "RS", "risk_group"), 
                     measure.vars = verified_clean_names,
                     variable.name = "Drug", value.name = "IC50")

# Lock the panel sequence to match your significance ranking exactly
scatter_data$Drug <- factor(scatter_data$Drug, levels = verified_clean_names)

# Generate a unified scatter plot grid for all 10 drugs (Full 2x5 Grid Layout)
p_scatter_10 <- ggscatter(scatter_data, x = "RS", y = "IC50",
                          add = "reg.line", conf.int = TRUE, 
                          add.params = list(color = "#B2182B", fill = "lightgray")) +
  stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top", size = 3) +
  facet_wrap(~Drug, scales = "free_y", ncol = 5) +
  labs(x = "Prognostic Risk Score (RS)", y = "Predicted ln(IC50)") +
  theme_bw() +
  theme(
    strip.background = element_rect(fill="gray95"),
    strip.text = element_text(face = "bold", size = 10)
  )

ggsave(filename = file.path(out_dir, "top10_drugs_IC50_vs_risk_score.pdf"), 
       plot = p_scatter_10, width = 14, height = 6.5)

