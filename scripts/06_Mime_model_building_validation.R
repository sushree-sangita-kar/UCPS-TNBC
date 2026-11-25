# 06  Mime: univariate Cox filter + 117 machine-learning combinations (mode = 'all'),
#     final RSF + CoxBoost model, C-index / KM / time-dependent AUC / ROC,
#     univariate Cox meta-analysis across cohorts, TCGA risk groups
#
# Input : results/objects/187tcga_list_training_validation_Datasets.Rdata (script 04c)
#         results/objects/candidates_malignant_only_187tcga.txt          (script 05)
# Output: results/objects/malignant_only_mime_result_tcga187.RData       (object: res)
#         results/objects/tcga_risk_malignant_only_tcga187.RData         (object: tcga_risk)
#         results/06_Mime/figures

library(xgboost)
library(survminer)
library(randomForestSRC)
library(Mime1)
source(here::here("config.R"))
out_dir <- step_dir("06_Mime")

load(file.path(OBJ_DIR, "187tcga_list_training_validation_Datasets.Rdata"))   # list_train_vali_Data (script 04c)
candidate_malignant_only <- readLines(file.path(OBJ_DIR, "candidates_malignant_only_187tcga.txt"))
dim(list_train_vali_Data$tcga) 

# Step 1: Build models
res <- ML.Dev.Prog.Sig(train_data = list_train_vali_Data$tcga,
                       list_train_vali_Data = list_train_vali_Data,
                       unicox.filter.for.candi = TRUE,
                       unicox_p_cutoff = 0.05,
                       candidate_genes = candidate_malignant_only,
                       mode = 'all',
                       nodesize = 5,
                       cores_for_parallel = 1,  
                       seed = 42)
  save(res, file = file.path(OBJ_DIR, "malignant_only_mime_result_tcga187.RData"))

# Step 2: C-index plot for all models (used to choose the final model)
p <- cindex_dis_all(res,
                    validate_set = names(list_train_vali_Data)[-1],
                    order = names(list_train_vali_Data),
                    width = 0.9)
ggsave(file.path(out_dir, "cindex_plot_tcga187.png"), plot = p,  width = 8, height = 10, dpi = 300)
ggsave(file.path(out_dir, "cindex_plot_tcga187.svg"), plot = p, width = 6, height = 10, dpi = 300)

# Step 3: final model chosen from the C-index plot
best_model <- "RSF + CoxBoost"

# Step 4: Survival curves — match loop to actual number of datasets
n_datasets <- length(names(list_train_vali_Data))
survplot <- vector("list", n_datasets)

for (i in 1:n_datasets) {
  survplot[[i]] <- rs_sur(res,
                          model_name = best_model,
                          dataset    = names(list_train_vali_Data)[i],
                          color      = c("blue", "red"),
                          median.line = "none",
                          cutoff     = 0.5,
                          conf.int   = TRUE,
                          xlab       = "Days",
                          pval.coord = c(100, 0.1)) +
    theme(
      plot.title= element_text(size = 28),
      axis.title.x = element_text(size = 24),   # "Days" label
      axis.title.y = element_text(size = 24),   # "Survival probability" label
      axis.text.x  = element_text(size = 22),   # x tick numbers (0, 2000, 4000…)
      axis.text.y  = element_text(size = 22),   # y tick numbers (0.00, 0.25…)
      legend.text  = element_text(size = 22),   # legend entries (Low n=94, High n=93)
      legend.title = element_text(size = 22)    # legend heading (Riskscore in tcga)
    )
  # enlarge the p / HR / CI annotation (a text layer inside the plot)
  for (j in seq_along(survplot[[i]]$layers)) {
    if (inherits(survplot[[i]]$layers[[j]]$geom, "GeomText")) {
      survplot[[i]]$layers[[j]]$aes_params$size <- 7
    }
  }
  print(survplot[[i]])
}

survival_plot <- aplot::plot_list(gglist = survplot, ncol = 3)
ggsave(file.path(out_dir, "KM_survival_plots_tcga187.svg"), plot = survival_plot,
       width = 25, height = 6, dpi = 300)

# Step 6: time-dependent AUC at 3 and 5 years (nearest-neighbour estimation, NNE)

# 3-year AUC
all.auc.3y <- cal_AUC_ml_res(res.by.ML.Dev.Prog.Sig = res,
                             train_data = list_train_vali_Data[["tcga"]],
                             inputmatrix.list = list_train_vali_Data,
                             mode = 'all',
                             AUC_time = 3,
                             auc_cal_method = "NNE")
# 5-year AUC
all.auc.5y <- cal_AUC_ml_res(res.by.ML.Dev.Prog.Sig = res,
                             train_data = list_train_vali_Data[["tcga"]],
                             inputmatrix.list = list_train_vali_Data,
                             mode = 'all',
                             AUC_time = 5,
                             auc_cal_method = "NNE")


# Step 7: AUC of all models
p_auc_3y <- auc_dis_all(all.auc.3y,
                        dataset = names(list_train_vali_Data),
                        validate_set = names(list_train_vali_Data)[-1],
                        order = names(list_train_vali_Data),
                        width = 0.35, year = 3)
ggsave(file.path(out_dir, "AUC_3year_tcga187.svg"), plot = p_auc_3y, width = 10, height = 12, dpi = 300)

p_auc_5y <- auc_dis_all(all.auc.5y,
                        dataset = names(list_train_vali_Data),
                        validate_set = names(list_train_vali_Data)[-1],
                        order = names(list_train_vali_Data),
                        width = 0.35, year = 5)
ggsave(file.path(out_dir, "AUC_5year_tcga187.svg"), plot = p_auc_5y, width = 10, height = 12, dpi = 300)

# Step 8: ROC curves of the final model
p_roc_3y <- roc_vis(all.auc.3y,
                    model_name = best_model,
                    dataset = names(list_train_vali_Data),
                    order = names(list_train_vali_Data),
                    anno_position = c(0.65, 0.55),
                    year = 3) +
  theme(axis.title = element_text(size = 26),
        axis.text = element_text(size = 30),
        legend.text = element_text(size = 32),
        legend.title = element_text(size = 32),
        strip.text = element_text(size = 26),
        plot.title = element_text(size = 28))

ggsave(file.path(out_dir, "ROC_3year.svg"), plot = p_roc_3y, width = 8.5, height = 8, dpi = 300)

p_roc_5y <- roc_vis(all.auc.5y,
                    model_name = best_model,
                    dataset = names(list_train_vali_Data),
                    order = names(list_train_vali_Data),
                    anno_position = c(0.65, 0.55),
                    year = 5) +
  theme(axis.title = element_text(size = 26),
        axis.text = element_text(size = 30),
        legend.text = element_text(size = 32),
        legend.title = element_text(size = 32),
        strip.text = element_text(size = 26),
        plot.title = element_text(size = 28))
ggsave(file.path(out_dir, "ROC_5year.svg"), plot = p_roc_5y, width = 8.5, height = 8, dpi = 300)


# Step 9: Unicox meta-analysis
unicox.rs.res <- cal_unicox_ml_res(res.by.ML.Dev.Prog.Sig = res,
                                   optimal.model = best_model,
                                   type = 'categorical')

metamodel <- cal_unicox_meta_ml_res(input = unicox.rs.res)

# Reorder the data frame columns to the order meta_unicox_vis expects
metamodel$data <- metamodel$data[, c(1:6, 8, 9, 7)]

# Generate the plot
p_meta <- meta_unicox_vis(metamodel,
                          dataset = names(list_train_vali_Data))


ggsave(file.path(out_dir, "meta_forest_plot.svg"), plot = p_meta, 
       width = 8.5, height = 2.5, dpi = 300, bg = "white")


# risk group stratification
str(res$riskscore, max.level = 2)
names(res$riskscore)
#Get risk scores for each dataset for the final model (RSF + CoxBoost)
risk_data <- res$riskscore$`RSF + CoxBoost`
# View the structure
str(risk_data)
#TCGA
tcga_risk <- risk_data$tcga
tcga_risk$risk_group <- ifelse(tcga_risk$RS > median(tcga_risk$RS), "High", "Low")
table(tcga_risk$risk_group)
save(tcga_risk, file = file.path(OBJ_DIR, "tcga_risk_malignant_only_tcga187.RData"))
mod <- res$ml.res[["RSF + CoxBoost"]]
class(mod)
# If mod is a list wrapping the fit, find the CoxBoost element first:
str(mod, max.level = 1)
## 1. Get exact names and coefficients
beta        <- as.numeric(coef(mod))
names(beta) <- mod$xnames

# Ensure these match print(mod$xnames) exactly
genes6      <- c("HELLS"  , "PUF60"  , "ASPHD1" , "KRT16"  , "PSMD4"  , "EXOSC4"  ,"RARRES1")
coef_std    <- beta[genes6]
## 2. Coefficients on the log2(TPM + 1) scale (stored coefficient / training SD)
train_sds   <- mod$sdx[match(genes6, mod$xnames)]
paper_betas <- coef_std / train_sds
names(paper_betas) <- genes6

print(paper_betas)


# Event counts per risk group per cohort
for (co in c("tcga","gse58812","gse96058")) {
  d <- res$riskscore[["RSF + CoxBoost"]][[co]]
  d$grp <- ifelse(d$RS > median(d$RS), "High","Low")
  print(co); print(table(group=d$grp, event=d$OS))
}

