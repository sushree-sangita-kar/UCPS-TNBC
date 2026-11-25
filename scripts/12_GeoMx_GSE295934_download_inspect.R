# 12  GSE295934 (GeoMx WTA, TNBC: >15 y survivors vs <3 y deaths)
#     Download the count matrix, probe annotation and sample metadata, and
#     check how count-matrix columns match the metadata.
#
# Output: data/GSE295934/ raw counts, PKC file, GEO series matrices,
#         GSE295934_meta.rds (sample metadata used by scripts 13 and 14)

suppressPackageStartupMessages({
  library(GEOquery)
  library(data.table)
  library(dplyr)
})

source(here::here("config.R"))
result_path <- GEOMX_DIR
dir.create(result_path, recursive = TRUE, showWarnings = FALSE)
options(timeout = 3600)

# 1. Supplementary files: gene-level raw counts + probe annotation (PKC)
base_url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE295nnn/GSE295934/suppl/"
supp_files <- c("GSE295934_raw_counts.tsv.gz",
                "GSE295934_Hs_R_NGS_WTA_v1.0.pkc.gz")

for (f in supp_files) {
  dest <- file.path(result_path, f)
  if (!file.exists(dest)) download.file(paste0(base_url, f), dest, mode = "wb")
}

counts <- fread(file.path(result_path, "GSE295934_raw_counts.tsv.gz"))

cat("\n========== RAW COUNTS ==========\n")
cat("Dimensions (rows x cols):", dim(counts), "\n")
cat("First column name:", colnames(counts)[1], "\n")
print(counts[1:5, 1:min(6, ncol(counts)), with = FALSE])
cat("\nFirst 10 column names:\n"); print(head(colnames(counts), 10))
cat("\nLast 5 column names:\n");   print(tail(colnames(counts), 5))
cat("\nRows that look like negative probes:\n")
print(grep("neg", counts[[1]], ignore.case = TRUE, value = TRUE))

# 2. Sample metadata from both SubSeries (GSE295858 + GSE295927)
meta <- bind_rows(lapply(c("GSE295858", "GSE295927"), function(g) {
  gse <- getGEO(g, GSEMatrix = TRUE, destdir = result_path)
  bind_rows(lapply(gse, function(x) pData(x)))
}))

cat("\n========== SAMPLE METADATA ==========\n")
cat("Number of samples in metadata:", nrow(meta), "\n")

# Compact summary of every metadata column: n unique values + first 6 values
cat("\nColumn summary (n_unique | first values):\n")
for (cc in colnames(meta)) {
  u <- unique(as.character(meta[[cc]]))
  cat(sprintf("- %s  [%d unique]: %s\n", cc, length(u),
              paste(head(u, 6), collapse = " || ")))
}

# Patients and segment types from the sample title (e.g. "TNBC001|ROI-019|PanCK-")
meta$patient <- sub("\\|.*$", "", meta$title)
meta$segment <- sub("^.*\\|", "", meta$title)
cat("\nSegments per patient x segment type (non-patient rows = likely NTC wells):\n")
print(table(meta$patient, meta$segment))

# 3. How do count-matrix columns match the metadata?
cnt_cols <- colnames(counts)[-1]
dcc_names <- basename(as.character(meta$supplementary_file_1))
cat("\n========== MATCHING ==========\n")
cat("Count columns matching GSM IDs      :", sum(cnt_cols %in% meta$geo_accession), "\n")
cat("Count columns matching sample titles:", sum(cnt_cols %in% meta$title), "\n")
cat("Count columns matching DCC names    :",
    sum(cnt_cols %in% dcc_names | cnt_cols %in% sub("\\.dcc(\\.gz)?$", "", dcc_names)), "\n")

saveRDS(meta, file.path(result_path, "GSE295934_meta.rds"))
cat("\nSaved metadata to:", file.path(result_path, "GSE295934_meta.rds"), "\n")