# 14  GSE295934 (GeoMx WTA, treatment-naive TNBC): analyses reported in the paper
#   17 long survivors (>15 y, GPx) vs 15 early deaths (<3 y, PPx)
#
#   (1) Locked 6-gene model in PanCK+ tumour epithelium
#       -> median risk by outcome, Wilcoxon test, AUC with DeLong 95% CI
#   (2) Tumour risk score vs composition of the adjacent stroma (SpatialDecon)
#       -> Spearman with stromal immune fraction; immune fraction High vs Low risk
#

# Standards followed
#  - Segments: the authors' QC-passed set.
#  - Q3 normalisation (GeomxTools RNA-NGS workflow; also used by the authors).
#    Q3 factors computed on the authors' filtered gene set, applied to all genes.
#  - Locked model, no refitting: predict(type = "lp") on patient-level PanCK+
#    pseudo-bulk, log2(Q3 + 1), the scale the model was trained on in Mime.
#  - Deconvolution: SpatialDecon (Danaher et al., Nat Commun 2022), set up as in
#    the Bioconductor NSCLC vignette: linear-scale Q3 data; per-segment background
#    from the negative probe (derive_GeoMx_background); raw counts for weights;
#    safeTME + safeTME.matches; PanCK+ segments given as pure tumour so cancer-cell
#    signal in stroma segments is not counted as TME cells.
#    Readout: prop_of_nontumor in PanCK- segments, averaged per patient.
#  - Immune fraction (pre-specified single summary): sum of all immune cell-type
#    proportions, i.e. every safeTME type except fibroblasts and endothelial cells.

suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(tidyr); library(pROC)
  library(CoxBoost); library(SpatialDecon); library(ggplot2); library(svglite)
})

# USER SETTINGS
source(here::here("config.R"))
result_path     <- GEOMX_DIR
counts_file     <- file.path(result_path, "GSE295934_raw_counts_fullWTA.tsv.gz")  # from script 13
authors_file    <- file.path(result_path, "GSE295934_raw_counts.tsv.gz")          # filtered genes + NegProbe-WTX
meta_file       <- file.path(result_path, "GSE295934_meta.rds")
res_path        <- file.path(OBJ_DIR, "malignant_only_mime_result_tcga187.RData")   # script 06
res_object_name <- "res"
model_name      <- "RSF + CoxBoost"
out_path        <- step_dir("14_GeoMx_GSE295934")
dir.create(out_path, recursive = TRUE, showWarnings = FALSE)
set.seed(1)

pal_outcome <- c(GPx = "#7FB77E", PPx = "#B35C9f")   # outcome groups: teal / purple
pal_risk    <- c(Low = "#1F77B4", High = "#C0392B")  # risk groups: unchanged


save_plot <- function(p, name, w = 4, h = 4) {
  ggsave(file.path(out_path, paste0(name, ".svg")), p, width = w, height = h, device = svglite)
  ggsave(file.path(out_path, paste0(name, ".png")), p, width = w, height = h, dpi = 300)
}

# Q3 normalisation as GeomxTools normalize(norm_method = "quant", desiredQuantile = .75):
# factor = Q3 of sample / geometric mean of all Q3s; normalised = counts / factor.
q3_norm <- function(m, q3_genes) {
  q3 <- apply(m[rownames(m) %in% q3_genes, , drop = FALSE], 2, quantile, probs = 0.75)
  sweep(m, 2, q3 / exp(mean(log(q3))), "/")
}

# 1. Load counts, negative probe and metadata
counts <- fread(counts_file)
cm <- as.matrix(counts[, -1, with = FALSE]); rownames(cm) <- counts[[1]]
auth <- fread(authors_file)
am <- as.matrix(auth[, -1, with = FALSE]); rownames(am) <- auth[[1]]
if (!"NegProbe-WTX" %in% rownames(am))
  stop("NegProbe-WTX row not found in the authors' matrix, so background cannot be derived.")
q3_genes <- setdiff(rownames(am), "NegProbe-WTX")

meta <- readRDS(meta_file) %>%
  mutate(library   = sub("^Library name:\\s*", "", description),
         patient   = sub("\\|.*$", "", title),
         segment   = sub("^.*\\|", "", title),
         prognosis = sub("^prognosis:\\s*", "", characteristics_ch1.1)) %>%
  filter(grepl("^TNBC", patient), prognosis %in% c("GPx", "PPx"))

keep    <- Reduce(intersect, list(colnames(cm), colnames(am), meta$library))
meta    <- meta[match(keep, meta$library), ]
cm      <- cm[, keep]
neg_raw <- am["NegProbe-WTX", keep]
meta$segment   <- factor(meta$segment,   levels = c("PanCK-", "PanCK+"))
meta$prognosis <- factor(meta$prognosis, levels = c("GPx", "PPx"))

pt <- distinct(meta, patient, prognosis)
cat("Segments:", ncol(cm), "| genes (full WTA):", nrow(cm), "\n")
cat("Patients by outcome:\n"); print(table(pt$prognosis))
cat("Segments by type:\n");    print(table(meta$segment))

# (1) LOCKED MODEL IN PanCK+ TUMOUR EPITHELIUM
pos <- meta$segment == "PanCK+"
pb  <- sapply(split(meta$library[pos], as.character(meta$patient[pos])),
              function(ids) rowSums(cm[, ids, drop = FALSE]))
pb_log <- log2(q3_norm(pb, q3_genes) + 1)

if (grepl("\\.rds$", res_path, ignore.case = TRUE)) {
  res <- readRDS(res_path)
} else {
  e <- new.env(); load(res_path, envir = e); res <- get(res_object_name, envir = e)
}
fit <- res$ml.res[[model_name]]
model_genes <- fit$xnames
beta_std    <- as.numeric(fit$coefficients[nrow(fit$coefficients), ])
nz_genes    <- model_genes[beta_std != 0]
cat("\nGenes with non-zero coefficient:", paste(nz_genes, collapse = ", "), "\n")

# Mime converts '-' to '.' in gene names, so harmonise before matching
X <- t(pb_log); colnames(X) <- gsub("-", ".", colnames(X))
missing_nz <- setdiff(nz_genes, colnames(X))
if (length(missing_nz) > 0)
  stop("Model genes absent from the matrix: ", paste(missing_nz, collapse = ", "),
       " -> run script 13 to rebuild the full WTA matrix.")
# genes the model carries with zero coefficient: placeholder at the training mean
for (g in setdiff(model_genes, colnames(X))) {
  X <- cbind(X, fit$meanx[match(g, model_genes)]); colnames(X)[ncol(X)] <- g
}
X <- X[, model_genes, drop = FALSE]

# predict() centres/scales with the TCGA training mean/SD, then applies the
# stored coefficients -> identical to the risk score used in TCGA
df <- data.frame(patient = rownames(X),
                 risk = as.numeric(predict(fit, newdata = X, type = "lp"))) %>%
  left_join(pt, by = "patient")
df$risk_group <- factor(ifelse(df$risk > median(df$risk), "High", "Low"), levels = c("Low", "High"))

wt      <- wilcox.test(risk ~ prognosis, data = df)
roc_obj <- roc(df$prognosis, df$risk, levels = c("GPx", "PPx"), direction = "<", quiet = TRUE)
auc_ci  <- ci.auc(roc_obj, method = "delong")
med_gpx <- median(df$risk[df$prognosis == "GPx"])
med_ppx <- median(df$risk[df$prognosis == "PPx"])

cat("\n===== (1) LOCKED MODEL IN PanCK+ TUMOUR EPITHELIUM =====\n")
cat(sprintf("Median risk  GPx = %.2f | PPx = %.2f | Wilcoxon p = %.3g\n", med_gpx, med_ppx, wt$p.value))
cat(sprintf("AUC = %.2f (95%% DeLong CI %.2f-%.2f)\n", auc_ci[2], auc_ci[1], auc_ci[3]))
fwrite(df, file.path(out_path, "patient_risk_scores_PanCKpos.csv"))

lab_outcome <- c(GPx = "Long-term survivors (GPx, >15 yr)", PPx = "Early deaths (PPx, <3 yr)")

p1 <- ggplot(df, aes(prognosis, risk, fill = prognosis)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.12, size = 1.6, show.legend = FALSE) +
  scale_fill_manual(values = pal_outcome, labels = lab_outcome, name = NULL) +
  labs(x = NULL, y = "Risk score (PanCK+ epithelial cells)",
       subtitle = sprintf("Wilcoxon p = %.3g", wt$p.value)) +
  theme_classic(base_size = 11) + theme(legend.position = "right")
save_plot(p1, "Fig_risk_score_by_outcome", w = 5.5)

p2 <- ggroc(roc_obj, legacy.axes = TRUE) +
  geom_abline(linetype = "dashed", colour = "grey60") +
  labs(subtitle = sprintf("AUC = %.2f (95%% CI %.2f-%.2f)", auc_ci[2], auc_ci[1], auc_ci[3])) +
  theme_classic(base_size = 11)
save_plot(p2, "Fig_ROC")

# (2) TUMOUR RISK SCORE vs ADJACENT STROMA COMPOSITION
# Normalise (Q3, linear scale) with the negative probe scaled by the same factors
genes    <- intersect(q3_genes, rownames(cm))
m_all    <- rbind(cm[genes, ], `NegProbe-WTX` = neg_raw)
norm_all <- q3_norm(m_all, q3_genes)
bg_all   <- derive_GeoMx_background(norm = norm_all, probepool = rep(1, nrow(norm_all)),
                                    negnames = "NegProbe-WTX")     # WTA = one probe pool
norm <- norm_all[genes, ]; bg <- bg_all[genes, ]; raw <- cm[genes, ]

data("safeTME"); data("safeTME.matches")
cat("\nsafeTME genes present in the data:", sum(rownames(safeTME) %in% genes),
    "of", nrow(safeTME), "\n")

# Deconvolve all segments together; PanCK+ segments = pure tumour
decon <- spatialdecon(norm = norm, bg = bg, X = safeTME, raw = raw,
                      cellmerges = safeTME.matches,
                      is_pure_tumor = meta$segment == "PanCK+",
                      n_tumor_clusters = 10)
saveRDS(decon, file.path(out_path, "spatialdecon_result.rds"))

orient <- function(m) { m <- as.matrix(m); if (nrow(m) == ncol(norm)) m else t(m) }  # segments x cells
prop <- orient(decon$prop_of_nontumor)
prop <- prop[, !grepl("^tumor", colnames(prop)), drop = FALSE]
prop[is.na(prop)] <- 0
rownames(prop) <- colnames(norm)
cell_types <- colnames(prop)

# Patient-level stromal composition: mean over that patient's PanCK- segments
neg <- meta$segment == "PanCK-"
pt_prop <- as.data.frame(prop[neg, , drop = FALSE]) %>%
  mutate(patient = as.character(meta$patient[neg])) %>%
  group_by(patient) %>%
  summarise(across(everything(), mean), n_stroma_segments = n(), .groups = "drop")
immune_types <- setdiff(cell_types, c("fibroblasts", "endothelial.cells"))
pt_prop$immune_fraction <- rowSums(pt_prop[, immune_types])
cat("Immune cell types summed:", paste(immune_types, collapse = ", "), "\n")

d <- inner_join(df[, c("patient", "risk", "risk_group", "prognosis")], pt_prop, by = "patient")
fwrite(d, file.path(out_path, "patient_risk_and_stroma_composition.csv"))

ct <- cor.test(d$risk, d$immune_fraction, method = "spearman", exact = FALSE)
wg <- wilcox.test(immune_fraction ~ risk_group, data = d, exact = FALSE)
med_low  <- median(d$immune_fraction[d$risk_group == "Low"])
med_high <- median(d$immune_fraction[d$risk_group == "High"])

cat("\n===== (2) TUMOUR RISK SCORE vs STROMAL IMMUNE FRACTION =====\n")
cat("Patients:", nrow(d), "\n")
cat(sprintf("Spearman rho = %.2f | p = %.3g\n", ct$estimate, ct$p.value))
cat(sprintf("Immune fraction  Low risk = %.2f | High risk = %.2f | Wilcoxon p = %.3g\n",
            med_low, med_high, wg$p.value))

p3 <- ggplot(d, aes(risk, immune_fraction)) +
  geom_point(aes(colour = prognosis), size = 2) +
  scale_colour_manual(values = pal_outcome) +
  labs(x = "Tumour (PanCK+) risk score", y = "Stromal immune fraction", colour = NULL,
       subtitle = sprintf("Spearman rho = %.2f, p = %.2g", ct$estimate, ct$p.value)) +
  theme_classic(base_size = 11)
save_plot(p3, "Fig_risk_vs_stromal_immune_fraction", w = 4.8)


p4 <- ggplot(d, aes(risk_group, immune_fraction, fill = risk_group)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) + geom_jitter(width = 0.12, size = 1.6) +
  scale_fill_manual(values = pal_risk) +
  scale_x_discrete(labels = c(Low = "Low risk", High = "High risk")) +
  labs(x = "Tumour risk group (median split)", y = "Stromal immune fraction",
       subtitle = sprintf("Wilcoxon p = %.2g", wg$p.value)) +
  theme_classic(base_size = 11) + theme(legend.position = "none")
save_plot(p4, "Fig_stromal_immune_fraction_by_risk_group")


