# 13  GSE295934: rebuild the full WTA gene-count matrix from DCC files
#
# Why: the authors' raw_counts.tsv is gene-filtered (10,150 genes) and drops
# ASPHD1, a non-zero gene in the locked model. The DCC files on GEO hold counts
# for every WTA probe, so the full matrix can be rebuilt without imputation.
#
# What it does
#  1. Downloads GSE295934_RAW.tar and extracts only the .dcc files
#  2. Reads probe annotation from the PKC (GeomxTools::readPKCFile)
#  3. Reads <Code_Summary> counts from every DCC
#  4. WTA has one probe per endogenous gene (GeomxTools vignette), so the gene
#     count = that probe's count; negative probes are summarised separately
#  5. Per-segment LOQ = geomean(neg) x geoSD(neg)^2, minimum 2 (vignette default)
#  6. Keeps only the segments that passed the authors' QC
#  7. SANITY CHECK: rebuilt counts must reproduce the authors' matrix for the
#     genes they kept. If not, stop and report.

suppressPackageStartupMessages({
  library(data.table)
  if (!requireNamespace("GeomxTools", quietly = TRUE)) BiocManager::install("GeomxTools", update = FALSE)
  library(GeomxTools)
})

source(here::here("config.R"))
result_path  <- GEOMX_DIR
dcc_dir      <- file.path(result_path, "dcc")
tar_file     <- file.path(result_path, "GSE295934_RAW.tar")
authors_file <- file.path(result_path, "GSE295934_raw_counts.tsv.gz")
pkc_gz       <- file.path(result_path, "GSE295934_Hs_R_NGS_WTA_v1.0.pkc.gz")
dir.create(dcc_dir, recursive = TRUE, showWarnings = FALSE)
options(timeout = 7200)


# 1. Download + extract DCC files (~736 MB tar; PNGs are skipped)

if (!file.exists(tar_file)) {
  download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE295nnn/GSE295934/suppl/GSE295934_RAW.tar",
                tar_file, mode = "wb")
}
dcc_in_tar <- grep("\\.dcc(\\.gz)?$", untar(tar_file, list = TRUE), value = TRUE)
cat("DCC files in tar:", length(dcc_in_tar), "\n")
untar(tar_file, files = dcc_in_tar, exdir = dcc_dir)
dcc_files <- list.files(dcc_dir, pattern = "\\.dcc(\\.gz)?$", full.names = TRUE)

# Library name per DCC file: GEO file names are "GSM<id>_<library>.dcc.gz", and the
# authors' count-matrix columns are the <library> part (all 1,151 match)
lib_ids <- sub("^GSM[0-9]+_", "", sub("\\.dcc(\\.gz)?$", "", basename(dcc_files)))


# 2. Probe annotation

pkc_file <- sub("\\.gz$", "", pkc_gz)
if (!file.exists(pkc_file)) writeLines(readLines(gzfile(pkc_gz)), pkc_file)
pkc <- as.data.frame(readPKCFile(pkc_file))
cat("PKC probes:", nrow(pkc), "\n"); print(head(pkc[, c("RTS_ID", "Target", "CodeClass")]))
print(table(pkc$CodeClass))


# 3. Read DCC <Code_Summary> blocks (readLines opens .gz transparently)

read_dcc_counts <- function(f) {
  x <- readLines(f, warn = FALSE)
  s <- grep("<Code_Summary>", x, fixed = TRUE); e <- grep("</Code_Summary>", x, fixed = TRUE)
  empty <- setNames(numeric(0), character(0))
  if (length(s) != 1 || length(e) != 1 || e - s < 2) return(empty)
  cs <- trimws(x[(s + 1):(e - 1)])
  cs <- cs[grepl(",", cs, fixed = TRUE)]
  if (length(cs) == 0) return(empty)
  sp <- do.call(rbind, strsplit(cs, ",", fixed = TRUE))
  setNames(as.numeric(sp[, 2]), sp[, 1])
}
probe_mat <- matrix(0, nrow = nrow(pkc), ncol = length(dcc_files),
                    dimnames = list(pkc$RTS_ID, lib_ids))
empty_dcc <- character(0)
for (i in seq_along(dcc_files)) {
  cnt <- read_dcc_counts(dcc_files[i])
  if (length(cnt) == 0) { empty_dcc <- c(empty_dcc, lib_ids[i]); next }
  k <- intersect(names(cnt), rownames(probe_mat))
  probe_mat[k, i] <- cnt[k]
  if (i %% 200 == 0) cat("  read", i, "of", length(dcc_files), "DCCs\n")
}
cat("DCCs with an empty Code_Summary block:", length(empty_dcc), "\n"); print(empty_dcc)


# 4. Gene-level matrix + negative-probe background

is_neg  <- pkc$CodeClass == "Negative"
endo    <- !is_neg
stopifnot(!anyDuplicated(pkc$Target[endo]))     # one probe per endogenous gene in WTA
gene_mat <- probe_mat[endo, , drop = FALSE]
rownames(gene_mat) <- pkc$Target[endo]

lneg <- log(pmax(probe_mat[is_neg, , drop = FALSE], 1))
neg_geomean <- exp(colMeans(lneg))
neg_geosd   <- exp(apply(lneg, 2, sd))
loq <- pmax(2, neg_geomean * neg_geosd^2)        # cutoff = 2 geoSD, minLOQ = 2


# 5. Keep authors' QC-passed segments + sanity check against their matrix

auth <- fread(authors_file)
am <- as.matrix(auth[, -1, with = FALSE]); rownames(am) <- auth[[1]]
am <- am[rownames(am) != "NegProbe-WTX", ]

seg_keep <- intersect(colnames(am), colnames(gene_mat))
cat("\nAuthors' segments:", ncol(am), "| found in DCCs:", length(seg_keep), "\n")
gene_mat <- gene_mat[, seg_keep]; loq <- loq[seg_keep]

cg <- intersect(rownames(am), rownames(gene_mat))
d  <- am[cg, seg_keep] - gene_mat[cg, seg_keep]
cat("Genes compared:", length(cg), "| fraction of identical values:", round(mean(d == 0), 4), "\n")

# GeomxTools shiftCountsOne() turns zeros into 1; apply it if that is the only difference
if (any(d != 0)) {
  only_shift <- all(d[d != 0] == 1 & gene_mat[cg, seg_keep][d != 0] == 0)
  if (only_shift) {
    cat("Differences are only zero -> 1 shifts (GeomxTools shiftCountsOne); applying the same.\n")
    gene_mat[gene_mat == 0] <- 1
    d <- am[cg, seg_keep] - gene_mat[cg, seg_keep]
  }
}
if (any(d != 0)) {
  print(summary(as.vector(d[d != 0])))
  stop("Rebuilt counts do not reproduce the authors' matrix.")
}
cat("Sanity check PASSED: rebuilt counts reproduce the authors' matrix exactly for", length(cg), "genes.\n")

cat("\nFull WTA matrix:", nrow(gene_mat), "genes x", ncol(gene_mat), "segments\n")
for (g in c("ASPHD1", "HELLS", "PUF60", "KRT16", "EXOSC4", "RARRES1")) {
  cat(sprintf("  %-8s present: %-5s | in authors' filtered set: %s | segments > LOQ: %.1f%%\n",
              g, g %in% rownames(gene_mat), g %in% rownames(am),
              if (g %in% rownames(gene_mat)) 100 * mean(gene_mat[g, ] > loq) else NA))
}

out_counts <- data.table(V1 = rownames(gene_mat), gene_mat)
fwrite(out_counts, file.path(result_path, "GSE295934_raw_counts_fullWTA.tsv.gz"), sep = "\t")


loq <- pmax(neg_geomean * neg_geosd^2, 2)[seg_keep]
summary(loq)
for (g in c("ASPHD1", "HELLS", "PUF60", "KRT16", "EXOSC4", "RARRES1"))
  cat(sprintf("  %-8s segments > LOQ: %.1f%%\n", g, 100 * mean(gene_mat[g, ] > loq)))
saveRDS(loq, file.path(result_path, "GSE295934_LOQ_per_segment.rds"))


cat("\nSaved full matrix + LOQ to:", result_path, "\n")
