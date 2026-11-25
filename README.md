# UCPS-TNBC

The UPS-stratified cell-derived prognostic signature (UCPS) is a 6-gene prognostic model for triple-negative breast cancer (KRT16, ASPHD1, EXOSC4, PUF60, RARRES1, HELLS). It was built from single-cell RNA-seq (GSE176078) and TCGA bulk RNA-seq, and validated in independent bulk and GeoMx spatial cohorts.

An online risk calculator for the locked model is hosted on the CG&T Lab website
https://project.iith.ac.in/cgntlab/TNBC_risk_calculator/ (source: `web/tnbc_risk_index.php`).

---

## Repository structure

```
UCPS-TNBC/
├── README.md
├── LICENSE
├── config.R        shared paths, sourced by every script
├── .here           marks the project root for the here package
├── resources/      small inputs kept in the repository (see "Stored inputs")
│   ├── ups_genes.txt                          UPS gene list used for AUCell
│   ├── UPS_AUCell_scores_epithelial_cells.csv per-cell UPS scores used in the study
│   ├── inferCNV_cell_annotation.csv           cell labels used for inferCNV
│   └── malignant_classification_anchored.csv  malignant calls used in the study
├── scripts/        analysis scripts, run in numerical order
├── web/tnbc_risk_index.php   risk calculator web page (locked 6-gene model)
├── data/           public input data (not tracked; see "Data")
└── results/        all outputs (not tracked; created when scripts run)
    └── objects/    .RData / .rds files passed from one script to the next
```

## Data

All data are public.

| Dataset | Source | Place in |
|---|---|---|
| GSE176078 scRNA-seq (CID44041, CID4465, CID4495, CID44971, CID44991, CID4515) | GEO | `data/GSE176078/<sample>/` |
| TCGA-BRCA expression (STAR counts, TPM) and clinical data (BCR Biotab) | GDC, downloaded by script 04a | `results/objects/` |
| GSE58812 | GEO, downloaded by script 04b | `data/GSE58812/` |
| GSE96058 | GEO: series records downloaded by script 04b; supplementary file `GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz` downloaded manually | `data/GSE96058/` |
| GDSC2 training data (expression and drug response) | https://osf.io/c6tfx/ | `data/oncoPredict/Training Data/` |
| GSE295934 GeoMx WTA (17 long survivors, 15 early deaths) | GEO, downloaded by script 12 | `data/GSE295934/` |
| `hg38_gencode_v27.txt` (inferCNV gene-order file) | inferCNV / Trinity CTAT | `data/reference/` |


## Scripts

Run the scripts in numerical order. Objects passed between scripts are saved in `results/objects/`.

### Single-cell discovery (GSE176078)

| Script | Description |
|---|---|
| `01_scRNA_QC_clustering_annotation.R` | Quality control, Harmony integration, clustering, cell-type annotation and epithelial sub-clustering |
| `02_AUCell_UPS_scoring.R` | How the per-cell AUCell scores of the UPS gene set were computed (later scripts use the stored scores in `resources/`) |
| `03_inferCNV_run.R` | inferCNV per sample with the stored cell labels and non-epithelial cells as reference; malignant calls with the reference-anchored threshold (later scripts use the stored calls in `resources/`) |

### Bulk cohorts (TCGA and GEO)

| Script | Description |
|---|---|
| `04a_TCGA_download_logTPM.R` | Downloads TCGA-BRCA expression and clinical data; builds the gene-level log2(TPM + 1) matrix and the TNBC survival table |
| `04b_GEO_validation_download.R` | Downloads and processes GSE58812 and the TNBC subset of GSE96058 |
| `04c_build_training_validation_cohorts.R` | Builds the TCGA training cohort (n = 187) and the GSE58812 and GSE96058 validation cohorts |
| `04d_TCGA_bulk_DEG_limma.R` | limma-voom differential expression, 187 TNBC tumours vs solid tissue normal |

### Signature development

| Script | Description |
|---|---|
| `05_malignant_call_UPS_split_scDEG.R` | Median UPS split within malignant cells, UPS-High vs UPS-Low DEGs, intersection with bulk DEGs (candidate genes) |
| `06_Mime_model_building_validation.R` | Univariate Cox filter, Mime (117 combinations), final RSF + CoxBoost model; C-index, Kaplan-Meier, time-dependent AUC/ROC, univariate Cox meta-analysis; TCGA risk groups |
| `07_multivariate_Cox_nomogram_calibration.R` | Multivariate Cox (age, stage, risk score), nomogram, 3- and 5-year calibration |

### Risk-group analyses (TCGA)

| Script | Description |
|---|---|
| `08_risk_group_DE_GSEA.R` | High- vs low-risk differential expression (limma-trend) and GSEA |
| `09_immune_ESTIMATE_CIBERSORT_checkpoints.R` | ESTIMATE, CIBERSORTx (LM22, run on the CIBERSORTx web portal) and immune-checkpoint genes |
| `10_oncoPredict_drug_response.R` | Drug-response prediction with oncoPredict (GDSC2) |

### Cell-cell communication (GSE176078)

| Script | Description |
|---|---|
| `11_CellChat_UPS_groups.R` | CellChat between UPS-High / UPS-Low malignant cells and all TME cell types |

### Spatial validation (GSE295934, GeoMx)

| Script | Description |
|---|---|
| `12_GeoMx_GSE295934_download_inspect.R` | Downloads GSE295934 counts, probe annotation and sample metadata |
| `13_GeoMx_GSE295934_rebuild_counts.R` | Rebuilds the full WTA count matrix from the DCC files |
| `14_GeoMx_GSE295934_risk_score_SpatialDecon.R` | Locked model in PanCK+ epithelium by outcome; tumour risk score vs stromal immune fraction (SpatialDecon) |

### Main outputs

| Script | Main output | Used by |
|---|---|---|
| 01 | `seurat_all_annotated.RData`, `epithelial_cells.RData` | 02, 03, 05, 11 |
| 04a | `TCGA_BRCA_data.RData`, `log_tpm_matrix_median.RData`, `clinical_biotab.RData`, `bulk_survival.RData` | 04c, 04d, 07 |
| 04b | `GSE58812_processed.RData`, `gse96058_tnbc_expr_log2_FPKM_plus_0.1.RData`, `gse96058_tnbc_survival.RData` | 04c |
| 04c | `187tcga_list_training_validation_Datasets.Rdata` | 04d, 06-10 |
| 04d | `bulk_degs_187.RData` | 05 |
| 05 | `candidates_malignant_only_187tcga.txt`, `seu_mal_UPS_groups.RData` | 06, 11 |
| 06 | `malignant_only_mime_result_tcga187.RData`, `tcga_risk_malignant_only_tcga187.RData` | 07-10, 14 |
| 12 | `GSE295934_meta.rds` (in `data/GSE295934/`) | 13, 14 |
| 13 | `GSE295934_raw_counts_fullWTA.tsv.gz` (in `data/GSE295934/`) | 14 |

## Software

R 4.3.3 for all scripts except script 06 (Mime), which was run with R 4.4.0.

| Package | Version | Scripts |
|---|---|---|
| Seurat | 5.4.0 | 01, 05, 11 |
| harmony | 1.2.4 | 01 |
| AUCell | 1.24.0 | 02 |
| infercnv | 1.18.1 | 03 |
| TCGAbiolinks | 2.30.4 | 04a |
| limma | 3.58.1 | 04d, 08 |
| Mime1 | 0.0.0.9000 | 06 |
| survival | 3.8.6 | 07 |
| survminer | 0.5.2 | 06, 07 |
| rms | 6.8.0 | 07 |
| clusterProfiler | 4.10.1 | 08 |
| fgsea | 1.28.0 | 08 |
| msigdbr | 26.1.0 (MSigDB v2026.1.Hs) | 08 |
| estimate | 1.0.13 | 09 |
| oncoPredict | 1.3.1 | 10 |
| CellChat | 2.2.0.9001 | 11 |
| GeomxTools | 3.5.0 | 13 |
| SpatialDecon | 1.12.3 | 14 |

CIBERSORTx was run on the web portal (https://cibersortx.stanford.edu).

## Licence

MIT (see `LICENSE`). Copyright (c) 2026 CG&T Lab, IIT Hyderabad.

