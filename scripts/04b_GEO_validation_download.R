# 04b  GEO validation cohorts: download and processing
#      GSE58812 (TNBC, Affymetrix HG-U133 Plus 2.0) and the TNBC subset of GSE96058
#
# Input : GEO, downloaded by getGEO() below
#         data/GSE96058/GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz
#           (GEO supplementary file of GSE96058; download it into data/GSE96058/ first, see README)
# Output: results/objects/GSE58812_processed.RData                    (object: GSE58812_processed)
#         results/objects/gse96058_tnbc_survival.RData                (object: gse96058_tnbc_survival)
#         results/objects/gse96058_tnbc_expr_log2_FPKM_plus_0.1.RData (object: gse96058_tnbc_expr_log2_FPKM_plus_0.1)

library(GEOquery)
library(hgu133plus2.db)
library(annotate)
library(dplyr)
library(data.table)
library(tidyverse)
source(here::here("config.R"))
out_dir <- step_dir("04b_GEO_validation")

# GSE58812
gse58812_dir <- file.path(DATA_DIR, "GSE58812")
dir.create(gse58812_dir, showWarnings = FALSE, recursive = TRUE)

gse <- getGEO("GSE58812", GSEMatrix = TRUE, destdir = gse58812_dir)
gse_data <- gse[[1]]
gse_expr <- exprs(gse_data)

# Probe -> gene symbol
gene_symbols <- getSYMBOL(rownames(gse_expr), "hgu133plus2.db")

# Keep probes with a symbol
valid_probes <- !is.na(gene_symbols)
gse_expr_filtered <- gse_expr[valid_probes, ]
gene_symbols_filtered <- gene_symbols[valid_probes]

probe_to_gene <- data.frame(
  probe_id = names(gene_symbols_filtered),
  gene_symbol = gene_symbols_filtered,
  stringsAsFactors = FALSE
)
probe_to_gene <- probe_to_gene[probe_to_gene$gene_symbol != "" &
                                 !is.na(probe_to_gene$gene_symbol), ]
gse_expr_filtered <- gse_expr_filtered[probe_to_gene$probe_id, ]

# One probe per gene: the probe with the highest mean expression
probe_to_gene$mean_expr <- rowMeans(gse_expr_filtered)
best_probes <- probe_to_gene %>%
  group_by(gene_symbol) %>%
  slice_max(mean_expr, n = 1, with_ties = FALSE) %>%
  ungroup()

gse_expr_genes <- gse_expr_filtered[best_probes$probe_id, ]
rownames(gse_expr_genes) <- best_probes$gene_symbol

# GEO values are MAS5-normalised on the linear scale (GEO data processing:
# MAS 5.0; range about 0.05 to 74,700), so log2-transform here
gse_expr_genes_log2 <- log2(gse_expr_genes + 1)

# Clinical data
sample_info <- pData(gse_data)
clinical_data <- data.frame(
  sample_id = rownames(sample_info),
  age_at_diagnosis = as.numeric(sample_info$`age at diag:ch1`),
  er_status = sample_info$`er-ihc:ch1`,
  pr_status = sample_info$`pr-ihc:ch1`,
  her2_status = sample_info$`her2-ihc:ch1`,
  death_event = as.numeric(sample_info$`death:ch1`),
  os_days = as.numeric(sample_info$`os (days):ch1`),
  os_years = as.numeric(sample_info$`os (days):ch1`) / 365.25,
  metastasis = sample_info$`meta:ch1`,
  mfs_days = as.numeric(sample_info$`mfs (days):ch1`),
  mfs_years = as.numeric(sample_info$`mfs (days):ch1`) / 365.25,
  stringsAsFactors = FALSE
)
cat("GSE58812:", nrow(gse_expr_genes_log2), "genes x", ncol(gse_expr_genes_log2), "samples\n")

GSE58812_processed <- list(
  expression = gse_expr_genes_log2,
  clinical = clinical_data
)
save(GSE58812_processed, file = file.path(OBJ_DIR, "GSE58812_processed.RData"))

# GSE96058: expression matrix, clinical data, TNBC filter
gse96058_dir <- file.path(DATA_DIR, "GSE96058")
expr_file <- file.path(gse96058_dir,
                       "GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz")
if (!file.exists(expr_file))
  stop("Download the GEO supplementary file ", basename(expr_file), " into ", gse96058_dir, " (see README)")

expr_raw <- fread(expr_file, data.table = FALSE, header = TRUE)
genes <- expr_raw[[1]]
expr_mat <- as.matrix(expr_raw[, -1])
rownames(expr_mat) <- genes

# Clinical data from the GSM records
g <- GEOquery::getGEO("GSE96058", GSEMatrix = FALSE, destdir = gse96058_dir)
gsm_list <- g@gsms

# Clinical info from one GSM record
extract_clinical_from_gsm <- function(gsm) {
  chars <- gsm@header$characteristics_ch1
  clinical <- list()
  clinical$geo_accession <- gsm@header$geo_accession
  clinical$title <- gsm@header$title
  for (char in chars) {
    # split on the first ": "
    parts <- strsplit(char, ": ", fixed = TRUE)[[1]]
    if (length(parts) >= 2) {
      key <- parts[1]
      value <- paste(parts[-1], collapse = ": ")
      key <- gsub(" ", "_", key)
      clinical[[key]] <- value
    }
  }
  return(clinical)
}

# All GSM records -> one data frame
process_gsm_list <- function(gsm_list) {
  clinical_list <- lapply(gsm_list, extract_clinical_from_gsm)
  all_cols <- unique(unlist(lapply(clinical_list, names)))
  clinical_df <- map_dfr(clinical_list, function(x) {
    missing_cols <- setdiff(all_cols, names(x))
    for (col in missing_cols) {
      x[[col]] <- NA
    }
    as_tibble(x)
  })
  return(clinical_df)
}

# TNBC: ER-, PR-, HER2-
identify_tnbc <- function(clinical_df) {
  clinical_df %>%
    filter(
      er_status == "0" &
        pgr_status == "0" &
        her2_status == "0"
    )
}

clinical_df <- process_gsm_list(gsm_list)
cat("GSE96058 total samples:", nrow(clinical_df), "\n")

clinical_df <- clinical_df %>%
  mutate(
    age_at_diagnosis = as.numeric(age_at_diagnosis),
    tumor_size = as.numeric(tumor_size),
    er_status = as.character(er_status),
    pgr_status = as.character(pgr_status),
    her2_status = as.character(her2_status),
    ki67_status = as.character(ki67_status),
    overall_survival_days = as.numeric(overall_survival_days),
    overall_survival_event = as.numeric(overall_survival_event),
    endocrine_treated = as.character(endocrine_treated),
    chemo_treated = as.character(chemo_treated)
  )

tnbc_df <- identify_tnbc(clinical_df)
cat("TNBC samples (ER = 0, PgR = 0, HER2 = 0):", nrow(tnbc_df), "\n")

cat("\nPAM50 subtypes in TNBC:\n");     print(table(tnbc_df$pam50_subtype, useNA = "ifany"))
cat("\nNHG grade in TNBC:\n");          print(table(tnbc_df$nhg, useNA = "ifany"))
cat("\nLymph node status in TNBC:\n");  print(table(tnbc_df$lymph_node_status, useNA = "ifany"))
cat("\nKi67 status in TNBC:\n");        print(table(tnbc_df$ki67_status, useNA = "ifany"))
cat("\nChemotherapy in TNBC:\n");       print(table(tnbc_df$chemo_treated, useNA = "ifany"))
cat("\nEndocrine treatment in TNBC:\n"); print(table(tnbc_df$endocrine_treated, useNA = "ifany"))
cat("\nAge in TNBC:\n");                print(summary(tnbc_df$age_at_diagnosis))
cat("\nTumour size in TNBC:\n");        print(summary(tnbc_df$tumor_size))
cat("\nSurvival events in TNBC:\n");    print(table(tnbc_df$overall_survival_event, useNA = "ifany"))

# TNBC survival table for samples present in the expression matrix
gse96058_tnbc_survival <- tnbc_df %>%
  filter(title %in% colnames(expr_mat)) %>%
  select(
    sample_id = title,
    time = overall_survival_days,
    status = overall_survival_event
  ) %>%
  mutate(
    time = as.numeric(time),
    status = as.numeric(status)
  ) %>%
  filter(!is.na(time) & !is.na(status))
save(gse96058_tnbc_survival, file = file.path(OBJ_DIR, "gse96058_tnbc_survival.RData"))

# TNBC expression. The GEO file is already log2(FPKM + 0.1) (GEO data processing
# description; observed range -3.32 to 23.05), so no further transformation.
expr_tnbc <- expr_mat[, gse96058_tnbc_survival$sample_id]
gse96058_tnbc_expr_log2_FPKM_plus_0.1 <- expr_tnbc
save(gse96058_tnbc_expr_log2_FPKM_plus_0.1,
     file = file.path(OBJ_DIR, "gse96058_tnbc_expr_log2_FPKM_plus_0.1.RData"))
cat("GSE96058 TNBC:", nrow(expr_tnbc), "genes x", ncol(expr_tnbc), "samples\n")

