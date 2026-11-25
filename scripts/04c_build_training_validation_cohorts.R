# 04c  Training and validation cohorts for Mime
#      TCGA (training, n = 187), GSE58812 and GSE96058 (validation)
#      Format per cohort: ID, OS.time (days), OS (1 = death), then one column per gene
#
# Input : results/objects/bulk_survival.RData           (script 04a)
#         results/objects/log_tpm_matrix_median.RData   (script 04a)
#         results/objects/GSE58812_processed.RData      (script 04b)
#         results/objects/gse96058_tnbc_expr_log2_FPKM_plus_0.1.RData,
#         results/objects/gse96058_tnbc_survival.RData  (script 04b)
# Output: results/objects/187tcga_list_training_validation_Datasets.Rdata
#         (object: list_train_vali_Data)

library(dplyr)
source(here::here("config.R"))

# TCGA (training)
load(file.path(OBJ_DIR, "bulk_survival.RData"))           # bulk_survival
load(file.path(OBJ_DIR, "log_tpm_matrix_median.RData"))   # log_tpm_matrix_median
tcga_expr <- log_tpm_matrix_median
tcga_clin <- bulk_survival

# Truncate sample barcodes to the 12-character patient ID
colnames(tcga_expr) <- substr(colnames(tcga_expr), 1, 12)
tcga_clin$barcode   <- substr(tcga_clin$barcode, 1, 12)

# One sample per patient: 194  tumour samples -> 187 patients
# (checked: for all 187 patients the column kept is a primary-tumour sample in bulk_survival)
dup_patients <- unique(tcga_clin$barcode[duplicated(tcga_clin$barcode)])
cat("Removing duplicates for:", length(dup_patients), "patients\n")
tcga_clin <- tcga_clin[!duplicated(tcga_clin$barcode), ]
tcga_expr <- tcga_expr[, !duplicated(colnames(tcga_expr))]

# Intersect and reorder to guarantee alignment
common    <- intersect(tcga_clin$barcode, colnames(tcga_expr))
tcga_clin <- tcga_clin[match(common, tcga_clin$barcode), ]
tcga_expr <- tcga_expr[, common]
stopifnot(all(tcga_clin$barcode == colnames(tcga_expr)))
cat("Final samples:", ncol(tcga_expr), "| clinical rows:", nrow(tcga_clin), "\n")

names(tcga_clin)[names(tcga_clin) == "barcode"] <- "ID"
tcga_clin <- tcga_clin %>%
  mutate(
    OS.time = as.numeric(time),
    OS      = ifelse(grepl("Dead|1", status), 1, 0),
    dataset = "TCGA"
  ) %>%
  filter(OS.time > 0, !is.na(OS)) %>%
  select(ID, OS.time, OS)

# Samples as rows, genes as columns
tcga_clin_aligned <- tcga_clin[match(colnames(tcga_expr), tcga_clin$ID), ]
tcga_merged <- cbind(
  tcga_clin_aligned[, c("ID", "OS.time", "OS")],
  t(tcga_expr[, tcga_clin_aligned$ID, drop = FALSE])
)
print(tcga_merged[1:5, 1:5])
stopifnot(nrow(tcga_merged) == 187)   # published training cohort

list_train_vali_Data <- list()
list_train_vali_Data[["tcga"]] <- tcga_merged

# GSE58812 (validation)
load(file.path(OBJ_DIR, "GSE58812_processed.RData"))      # GSE58812_processed
gse58812_expr <- GSE58812_processed$expression
gse58812_clin <- GSE58812_processed$clinical
cat("GSE58812 expression:", dim(gse58812_expr), "| clinical:", dim(gse58812_clin), "\n")

gse58812_clin_clean <- gse58812_clin %>%
  mutate(
    ID      = sample_id,
    OS.time = as.numeric(os_days),
    OS      = as.numeric(death_event)
  ) %>%
  filter(OS.time > 0, !is.na(OS.time), !is.na(OS)) %>%
  select(ID, OS.time, OS)
cat("GSE58812 samples:", nrow(gse58812_clin_clean), "| events:", sum(gse58812_clin_clean$OS), "\n")

common_samples_58812  <- intersect(colnames(gse58812_expr), gse58812_clin_clean$ID)
gse58812_clin_aligned <- gse58812_clin_clean[match(common_samples_58812, gse58812_clin_clean$ID), ]
gse58812_expr_aligned <- gse58812_expr[, common_samples_58812, drop = FALSE]

gse58812_merged <- cbind(
  gse58812_clin_aligned[, c("ID", "OS.time", "OS")],
  t(gse58812_expr_aligned)
)
cat("GSE58812 merged:", dim(gse58812_merged), "\n")

# GSE96058 (validation)
load(file.path(OBJ_DIR, "gse96058_tnbc_expr_log2_FPKM_plus_0.1.RData"))
load(file.path(OBJ_DIR, "gse96058_tnbc_survival.RData"))

gse96058_expr <- gse96058_tnbc_expr_log2_FPKM_plus_0.1
gse96058_clin <- gse96058_tnbc_survival
gse96058_clin_clean <- gse96058_clin %>%
  mutate(
    ID     = sample_id,
    time   = as.numeric(time),
    status = as.numeric(status)
  ) %>%
  filter(time > 0, !is.na(time), !is.na(status)) %>%
  select(ID, time, status)
cat("GSE96058 samples:", nrow(gse96058_clin_clean), "| events:", sum(gse96058_clin_clean$status), "\n")

common_samples_96058  <- intersect(colnames(gse96058_expr), gse96058_clin_clean$ID)
gse96058_clin_aligned <- gse96058_clin_clean[match(common_samples_96058, gse96058_clin_clean$ID), ]
gse96058_expr_aligned <- gse96058_expr[, common_samples_96058, drop = FALSE]

gse96058_merged <- cbind(
  gse96058_clin_aligned[, c("ID", "time", "status")],
  t(gse96058_expr_aligned)
) %>%
  rename(OS = status, OS.time = time)
cat("GSE96058 merged:", dim(gse96058_merged), "\n")

# Combine and save
list_train_vali_Data[["gse58812"]] <- gse58812_merged
list_train_vali_Data[["gse96058"]] <- gse96058_merged

for (ds in names(list_train_vali_Data)) {
  df <- list_train_vali_Data[[ds]]
  cat(sprintf("%s: %d samples, %d genes, %d events (%.1f%%)\n",
              ds, nrow(df), ncol(df) - 3, sum(df$OS), 100 * sum(df$OS) / nrow(df)))
}

save(list_train_vali_Data,
     file = file.path(OBJ_DIR, "187tcga_list_training_validation_Datasets.Rdata"))
