# =============================================================================
# config.R - shared paths for every script in this repository
#
# Every script starts with source(here::here("config.R")).
# here() finds the project root through the empty ".here" file in the top
# folder, so the scripts run from any working directory.
# Change DATA_DIR only if the raw data are kept somewhere else.
# =============================================================================

suppressPackageStartupMessages(library(here))

DATA_DIR    <- here("data")       # downloaded public data (not tracked by git)
RES_DIR     <- here("resources")  # small input tables kept in the repository
RESULTS_DIR <- here("results")    # all outputs (not tracked by git)
OBJ_DIR     <- file.path(RESULTS_DIR, "objects")   # .RData/.rds passed between scripts

# Input data folders (see README: "Data")
SC_DIR     <- file.path(DATA_DIR, "GSE176078")     # one sub-folder per sample
REF_DIR    <- file.path(DATA_DIR, "reference")     # inferCNV gene-order file
GEOMX_DIR  <- file.path(DATA_DIR, "GSE295934")     # GeoMx download (script 13)

# Folder for one script's outputs, created on first use
step_dir <- function(name) {
  d <- file.path(RESULTS_DIR, name)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

dir.create(OBJ_DIR, recursive = TRUE, showWarnings = FALSE)
