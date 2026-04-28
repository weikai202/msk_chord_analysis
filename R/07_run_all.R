# =============================================================================
# 07_run_all.R — Master script: run entire pipeline in order
# =============================================================================
# Usage:
#   1. Download MSK-CHORD 2024 from cBioPortal (see README)
#   2. Place unzipped files in data/raw/
#   3. Rscript R/07_run_all.R   (or source() from RStudio)
# =============================================================================

library(here)

# ── Pipeline steps ─────────────────────────────────────────────────────────────
steps <- list(
  list(script = "R/01_data_processing.R",   label = "Step 1: Data processing & imputation"),
  list(script = "R/02_descriptive.R",        label = "Step 2: Descriptive analysis & Table 1"),
  list(script = "R/03_chi_square_fdr.R",     label = "Step 3: Chi-square tests + BH-FDR"),
  list(script = "R/04_logistic_regression.R",label = "Step 4: Multivariable logistic regression"),
  list(script = "R/05_survival_analysis.R",  label = "Step 5: Survival analysis (KM + Cox PH)"),
  list(script = "R/06_replication_cohorts.R",label = "Step 6: Replication in CRC & Breast")
)

for (step in steps) {
  message("\n", paste(rep("=", 60), collapse = ""))
  message("▶ ", step$label)
  message(paste(rep("=", 60), collapse = ""))
  tryCatch(
    source(here::here(step$script)),
    error = function(e) {
      message("ERROR in ", step$script, ": ", e$message)
      message("Pipeline continuing with next step...")
    }
  )
}

message("\n", paste(rep("✓", 30), collapse = ""))
message("Pipeline complete! Outputs in:")
message("  Figures: ", here::here("output", "figures"))
message("  Tables:  ", here::here("output", "tables"))
