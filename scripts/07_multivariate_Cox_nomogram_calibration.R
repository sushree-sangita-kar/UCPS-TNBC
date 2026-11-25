## 07  TCGA (n = 187): multivariate Cox (age, stage, UCPS risk score),
##     nomogram, calibration at 3 and 5 years
##
## Input : results/objects/187tcga_list_training_validation_Datasets.Rdata  (script 04c)
##         results/objects/clinical_biotab.RData   (script 04a)
##         results/objects/malignant_only_mime_result_tcga187.RData   (script 06)
## Output: results/07_multivariate_Cox_nomogram/

library(dplyr)
source(here::here("config.R"))
load(file.path(OBJ_DIR, "187tcga_list_training_validation_Datasets.Rdata"))
load(file.path(OBJ_DIR, "clinical_biotab.RData"))
result_path <- step_dir("07_multivariate_Cox_nomogram")
# radiation
radiation <- clinical_biotab[["clinical_radiation_brca"]]
radiation<- radiation[-c(1, 2), ]
colnames(radiation)
# drug chemo
drug <- clinical_biotab[["clinical_drug_brca"]]
drug <-drug [-c(1, 2), ]
colnames(drug)
# Stage
# patient id in traing dataset
patient_ID_training <- list_train_vali_Data$tcga$ID
# Extract 12-character patient IDs from training data
patient_ids <- substr(patient_ID_training, 1, 12)


patient <- clinical_biotab[["clinical_patient_brca"]]
patient<- patient[-c(1, 2), ]
colnames(patient)

# Filter clinical tables for your TNBC patients
patient_filtered <- patient[patient$bcr_patient_barcode %in% patient_ids, ]
cat("Matched patients:", nrow(patient_filtered), "\n")

# Clean and prepare
patient_filtered$age <- as.numeric(patient_filtered$age_at_diagnosis)

## check the raw stage labels first, so nothing is silently dropped
print(table(patient_filtered$ajcc_pathologic_tumor_stage, useNA = "ifany"))

# Simplify stage
## added "^Stage III$" (plain Stage III was becoming NA)
patient_filtered$stage <- case_when(
  grepl("^Stage IIIA|^Stage IIIB|^Stage IIIC|^Stage III$", patient_filtered$ajcc_pathologic_tumor_stage) ~ "III",
  grepl("^Stage IV", patient_filtered$ajcc_pathologic_tumor_stage) ~ "IV",
  grepl("^Stage IIA|^Stage IIB|^Stage II$", patient_filtered$ajcc_pathologic_tumor_stage) ~ "II",
  grepl("^Stage IA|^Stage IB|^Stage I$", patient_filtered$ajcc_pathologic_tumor_stage) ~ "I",
  TRUE ~ NA_character_
)
print(table(patient_filtered$stage, useNA = "ifany"))


# Multivariate forest plot

# Multivariate: direct ggforest

library(survival)
library(survminer)
library(ggplot2)
library(gridExtra)

load(file.path(OBJ_DIR, "malignant_only_mime_result_tcga187.RData"))   # res

# Step 1: Merge risk score with clinical data
# Get risk scores from the best model
risk_scores <- res$riskscore$`RSF + CoxBoost`$tcga


risk_scores$ID <- substr(risk_scores$ID, 1, 12)

# Merge
nomo_data <- merge(patient_filtered, risk_scores, 
                   by.x = "bcr_patient_barcode", by.y = "ID")


# Remove NAs
nomo_data <- nomo_data[!is.na(nomo_data$stage) & 
                         #!is.na(nomo_data$N_stage) & 
                         !is.na(nomo_data$age) &
                         !is.na(nomo_data$OS.time) &
                         nomo_data$OS.time > 0, ]

# Convert to factors
nomo_data$stage <- factor(nomo_data$stage, levels = c("I", "II", "III", "IV"))
#nomo_data$N_stage <- factor(nomo_data$N_stage, levels = c("N0", "N1", "N2", "N3"))

cat("Patients for nomogram:", nrow(nomo_data), "\n")
cat("Events:", sum(nomo_data$OS), "\n")


# 1. Prepare variables as numeric ordinal
stage_map <- c("Stage I" = 1, "Stage II" = 2, "Stage III" = 3, "Stage IV" = 4,
               "I" = 1, "II" = 2, "III" = 3, "IV" = 4,
               "1" = 1, "2" = 2, "3" = 3, "4" = 4)

nomo_data$Stage<- as.numeric(stage_map[as.character(nomo_data$stage )])
#nomo_data$N_stage <- as.numeric(nomo_data$N_stage)
nomo_data$Risk_score <- as.numeric(nomo_data$RS)
nomo_data$Age <- as.numeric(nomo_data$age)

# 2. Multivariate model (all factors)
cox_multi <- coxph(Surv(OS.time, OS) ~ Age + Stage + Risk_score,
                   data = nomo_data)


# Multivariate forest plot via ggforest
p_multi <- ggforest(cox_multi,
                    data = nomo_data,
                    main = "Multivariate Cox regression analysis",
                    cpositions = c(0.02, 0.22, 0.4),
                    fontsize = 1.6,        # ggforest uses cex-style scaling, NOT pt size
                    refLabel = "reference",
                    noDigits = 3)

# Save multivariate
ggsave(file.path(result_path, "forest_multivariate_ggforest.svg"),
       p_multi,
       width = 28, height = 16,    # inches — adjust based on number of variables
       dpi = 300)

# Nomogram
# Style: Blue density waves + red dots + significance stars
# All multivariate factors included

library(survival)
library(rms)
library(regplot)

# 2. Fit cph model (rms) with ALL multivariate factors
# Convert OS.time from days to years
nomo_data$OS.time_yr <- nomo_data$OS.time / 365.25

dd <- datadist(nomo_data)
options(datadist = "dd")

fit_nomo <- cph(Surv(OS.time_yr, OS) ~ Age + Stage + Risk_score,
                data = nomo_data,
                x = TRUE, y = TRUE, surv = TRUE,
                time.inc = 5)

print(fit_nomo)

# 3. regplot nomogram (blue waves + red dots + stars)

## Example patient shown on the nomogram: alive beyond 5 years, with a linear
## predictor (= total points) closest to the cohort median.
lp_nomo    <- predict(fit_nomo)                     # linear predictor
candidates <- which(nomo_data$OS.time_yr > 5 & nomo_data$OS == 0)
best       <- candidates[which.min(abs(lp_nomo[candidates] - median(lp_nomo)))]
print(nomo_data[best, c("bcr_patient_barcode", "Age", "Stage", "Risk_score", "OS.time_yr", "OS")])

## The axis label "Pr( OS.time_yr > 3 )" was relabelled "3-year OS" in Inkscape.
regplot(fit_nomo,
        observation = nomo_data[best, ],
        title    = "Nomogram",
        points   = TRUE,
        plots    = c("density", "boxes"),
        failtime = c(3, 5),                # years, matching the model's time unit
        prfail   = FALSE,                  # survival probability, not failure probability
        droplines = TRUE,
        rank     = "sd",
        subticks = TRUE,
        showP    = TRUE,
        interval = "confidence",
        leftlabel = TRUE,
        color.palette = c("#4292C6", "#EF3B2C"),
        cex = 0.85)

dev.copy(svglite::svglite, file.path(result_path, "nomogram_regplot.svg"),
         width = 8, height = 5)
dev.off()



# Calibration curves

library(survival)
library(rms)
library(ggplot2)

# 1. Fit full model (time in days)
dd <- datadist(nomo_data)
options(datadist = "dd")

fit_full <- cph(Surv(OS.time, OS) ~ Age + Stage + Risk_score,
                data = nomo_data,
                x = TRUE, y = TRUE, surv = TRUE,
                time.inc = 1095)

# 2. Bootstrap validation -> C-index
set.seed(42)
val <- validate(fit_full, method = "boot", B = 200)
print(val)

corrected_Dxy <- val["Dxy", "index.corrected"]
corrected_C   <- 0.5 + (corrected_Dxy / 2)

## 95% CI from the standard error of the C-index (Hmisc::rcorr.cens).
## S.D. is the SE of Dxy, so SE(C) = S.D./2. Negative linear predictor because
## a higher lp means higher risk (shorter survival).
rc         <- rcorr.cens(-predict(fit_full), Surv(nomo_data$OS.time, nomo_data$OS))
C_apparent <- unname(rc["C Index"])
se_C       <- unname(rc["S.D."]) / 2

C_lower <- max(corrected_C - 1.96 * se_C, 0)
C_upper <- min(corrected_C + 1.96 * se_C, 1)

cat(sprintf("Apparent C-index: %.3f (95%% CI %.3f-%.3f)\n",
            C_apparent, C_apparent - 1.96 * se_C, C_apparent + 1.96 * se_C))
cat(sprintf("Bootstrap-corrected C-index (B = 200): %.3f (95%% CI %.3f-%.3f)\n",
            corrected_C, C_lower, C_upper))

c_index_text <- paste0("C-index:\n",
                       sprintf("%.3f", corrected_C),
                       " (95% CI: ",
                       sprintf("%.3f", C_lower), "-",
                       sprintf("%.3f", C_upper), ")")

# 3. Calibration curves (3, 5-year)
## cmethod = "KM" (the standard rms calibration plot with error bars).
## m = subjects per group; about n/4 gives ~4 groups here.
m_grp <- round(nrow(nomo_data) / 4)

# 3-year
fit_3yr <- cph(Surv(OS.time, OS) ~ Age + Stage + Risk_score,
               data = nomo_data,
               x = TRUE, y = TRUE, surv = TRUE,
               time.inc = 1095)
set.seed(42)
cal_3yr <- calibrate(fit_3yr, cmethod = "KM", method = "boot",
                     u = 1095, m = m_grp, B = 200)

# 5-year
fit_5yr <- cph(Surv(OS.time, OS) ~ Age + Stage + Risk_score,
               data = nomo_data,
               x = TRUE, y = TRUE, surv = TRUE,
               time.inc = 1825)
set.seed(42)
cal_5yr <- calibrate(fit_5yr, cmethod = "KM", method = "boot",
                     u = 1825, m = m_grp, B = 200)

## Reference check: rms's own plot (X = bootstrap-corrected, bars = 95% CI of KM).
## The ggplot below should match these.
svglite::svglite(file.path(result_path, "calibration_rms_check.svg"), width = 10, height = 5)
par(mfrow = c(1, 2))
plot(cal_3yr, xlim = c(0, 1), ylim = c(0, 1), main = "3-year")
abline(0, 1, lty = 2, col = "grey60")
plot(cal_5yr, xlim = c(0, 1), ylim = c(0, 1), main = "5-year")
abline(0, 1, lty = 2, col = "grey60")
dev.off()

## extraction for the KM method. No extra binning, since calibrate()
## already groups patients. std.err is the SE of -log(KM) (rms::groupkm),
## so the 95% CI is KM * exp(-/+ 1.96 * SE).
extract_cal_km <- function(cal_obj, label) {
  needed <- c("mean.predicted", "KM", "KM.corrected", "std.err")
  if (!all(needed %in% colnames(cal_obj))) {
    stop("Unexpected calibrate() columns: ", paste(colnames(cal_obj), collapse = ", "))
  }
  z  <- qnorm(0.975)
  km <- cal_obj[, "KM"]
  se <- cal_obj[, "std.err"]
  data.frame(
    predicted = cal_obj[, "mean.predicted"],
    observed  = km,                          # observed KM (apparent)
    corrected = cal_obj[, "KM.corrected"],   # bootstrap-corrected
    lower     = pmax(km * exp(-z * se), 0),
    upper     = pmin(km * exp( z * se), 1),
    timepoint = label
  )
}

cal_df <- rbind(
  extract_cal_km(cal_3yr, "3-year"),
  extract_cal_km(cal_5yr, "5-year")
)
cal_df$timepoint <- factor(cal_df$timepoint, levels = c("3-year", "5-year"))
print(cal_df)

cal_colors <- c("3-year" = "#666666",
                "5-year" = "#56B4E9")

## same layout as rms's plot. Line, dots and bars = observed KM with
## 95% CI; X = bootstrap-corrected estimate. Styling kept from your version.
p_cal <- ggplot(cal_df, aes(x = predicted, color = timepoint)) +
  geom_abline(intercept = 0, slope = 1, linetype = "solid",
              color = "grey70", linewidth = 0.6) +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                width = 0.02, linewidth = 0.5) +
  geom_line(aes(y = observed), linewidth = 0.8) +
  geom_point(aes(y = observed), size = 2) +
  geom_point(aes(y = corrected), shape = 4, size = 3, stroke = 1.2) +
  scale_color_manual(values = cal_colors) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2),
                     name = "Nomogram-predicted OS probability") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2),
                     name = "Observed OS probability") +
  annotate("text",
           x = 0.55, y = 0.15,
           label = c_index_text,
           size = 4.5, fontface = "bold", hjust = 0) +
  theme_bw(base_size = 16) +
  theme(
    legend.position = c(0.15, 0.85),
    legend.background = element_rect(fill = "white", color = "grey80"),
    legend.title = element_blank(),
    legend.text = element_text(size = 16),
    axis.text = element_text(color = "black", size = 16),
    axis.title = element_text(size = 17),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey92")
  )

ggsave(file.path(result_path, "calibration_combined.svg"), p_cal,
       width = 7, height = 6, dpi = 300)

