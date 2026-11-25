# 04a  TCGA-BRCA download (TCGAbiolinks) and gene-level log2(TPM + 1) matrix
#
#
# Output: results/objects/TCGA_BRCA_data.RData         (object: data, SummarizedExperiment)
#         results/objects/clinical_biotab.RData        (object: clinical_biotab)
#         results/objects/log_tpm_matrix_median.RData  (object: log_tpm_matrix_median)
#         results/objects/bulk_survival.RData          (object: bulk_survival; TNBC
#                                                       female tumours with valid survival)

library(TCGAbiolinks)
library(SummarizedExperiment)
source(here::here("config.R"))
set.seed(1234)

## 1. Expression (STAR counts / TPM)
query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification"
)
GDCdownload(query = query)
data <- GDCprepare(query = query)
save(data, file = file.path(OBJ_DIR, "TCGA_BRCA_data.RData"))

## 2. Clinical (BCR Biotab)
query_clinical <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Clinical",
  data.type = "Clinical Supplement",
  data.format = "BCR Biotab"
)
GDCdownload(query_clinical)
clinical_biotab <- GDCprepare(query_clinical)
save(clinical_biotab, file = file.path(OBJ_DIR, "clinical_biotab.RData"))

## 3. log2(TPM + 1), gene symbols, median over duplicated symbols
tpm_matrix     <- assay(data, "tpm_unstrand")
log_tpm_matrix <- log2(tpm_matrix + 1)

gene_info <- rowData(data)
id_map <- data.frame(
  ensembl_id = gene_info$gene_id,
  gene_name  = gene_info$gene_name,
  stringsAsFactors = FALSE
)

# Remove entries where one row maps to multiple genes (gene_name contains ";")
single_gene_entries <- !grepl(";", id_map$gene_name)
id_map         <- id_map[single_gene_entries, ]
log_tpm_matrix <- log_tpm_matrix[single_gene_entries, ]
rownames(log_tpm_matrix) <- id_map$gene_name

# Median for gene symbols with more than one Ensembl ID
log_tpm_matrix_median <- aggregate(log_tpm_matrix,
                                   by = list(Gene = rownames(log_tpm_matrix)),
                                   FUN = median)
rownames(log_tpm_matrix_median) <- log_tpm_matrix_median$Gene
log_tpm_matrix_median <- as.matrix(log_tpm_matrix_median[, -1])

dim(log_tpm_matrix)         
dim(log_tpm_matrix_median) 
save(log_tpm_matrix_median, file = file.path(OBJ_DIR, "log_tpm_matrix_median.RData"))

## 4.female tumours with valid overall survival
data_tnbc <- data

colData(data_tnbc)$paper_BRCA_Subtype_PAM50[colData(data_tnbc)$paper_BRCA_Subtype_PAM50 == "Normal" &
                                              colData(data_tnbc)$tissue_type == "Tumor"] <- "Normal_like"
colData(data_tnbc)$paper_BRCA_Subtype_PAM50[is.na(colData(data_tnbc)$paper_BRCA_Subtype_PAM50) &
                                              colData(data_tnbc)$tissue_type == "Normal"] <- "Normal"
table(colData(data_tnbc)$paper_BRCA_Subtype_PAM50, useNA = "always")


data_tnbc <- data_tnbc[, !is.na(colData(data_tnbc)$gender) &
                         colData(data_tnbc)$gender == "female" &
                         !is.na(colData(data_tnbc)$paper_BRCA_Subtype_PAM50) &
                         colData(data_tnbc)$paper_BRCA_Subtype_PAM50 == "Basal"]

# Valid survival: Dead with days_to_death > 0, or Alive with days_to_last_follow_up > 0
data_tnbc <- data_tnbc[, !is.na(colData(data_tnbc)$vital_status) &
                         ((colData(data_tnbc)$vital_status == "Dead" & !is.na(colData(data_tnbc)$days_to_death) &
                             colData(data_tnbc)$days_to_death > 0) |
                            (colData(data_tnbc)$vital_status == "Alive" & !is.na(colData(data_tnbc)$days_to_last_follow_up) &
                               colData(data_tnbc)$days_to_last_follow_up > 0))]
table(colData(data_tnbc)$vital_status, useNA = "always")

sample_info <- as.data.frame(colData(data_tnbc))
bulk_survival <- data.frame(
  barcode = sample_info$barcode,
  time    = ifelse(!is.na(sample_info$days_to_death), sample_info$days_to_death,
                   sample_info$days_to_last_follow_up),
  status  = ifelse(!is.na(sample_info$days_to_death), 1, 0)   # 1 = death
)

save(bulk_survival, file = file.path(OBJ_DIR, "bulk_survival.RData"))

