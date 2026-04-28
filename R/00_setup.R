# =============================================================================
# 00_setup.R — Package installation and project configuration
# MSK-CHORD 2024 Racial/Ethnic Disparities in Driver Gene Mutations
# =============================================================================
# STATISTICAL RATIONALE (统计原理):
#   This script installs all required packages for reproducible analysis.
#   Setting a global random seed (set.seed(42)) ensures that stochastic
#   procedures — MICE imputation chains, bootstrap CIs — produce identical
#   results across machines and R sessions.
# =============================================================================

# ── 0.1 Install missing packages ─────────────────────────────────────────────
required_pkgs <- c(
  "tidyverse",   # data wrangling + ggplot2
  "broom",       # tidy model outputs
  "survival",    # Kaplan-Meier + Cox PH
  "survminer",   # publication-ready KM plots
  "mice",        # multiple imputation
  "patchwork",   # compose multi-panel plots
  "scales",      # axis formatting
  "gtsummary",   # Table 1 generation
  "rms",         # cox.zph Schoenfeld residuals
  "viridis",     # colour-blind-safe palettes
  "tableone",    # Table 1 alternative
  "flextable",   # export tables to Word/PDF
  "writexl",     # write Excel tables
  "here"         # portable file paths
)

new_pkgs <- required_pkgs[!required_pkgs %in% installed.packages()[, "Package"]]
if (length(new_pkgs) > 0) {
  install.packages(new_pkgs, repos = "https://cloud.r-project.org")
}

# Load all packages
invisible(lapply(required_pkgs, library, character.only = TRUE))

# ── 0.2 Global settings ───────────────────────────────────────────────────────
set.seed(42)          # reproducibility seed
options(scipen = 999) # suppress scientific notation in outputs

# ── 0.3 Directory paths ───────────────────────────────────────────────────────
# Adjust DATA_DIR to point to the unzipped MSK-CHORD 2024 folder
DATA_DIR   <- here::here("data", "raw")
OUTPUT_FIG <- here::here("output", "figures")
OUTPUT_TAB <- here::here("output", "tables")

# ── 0.4 Driver genes of interest ─────────────────────────────────────────────
DRIVER_GENES <- c("EGFR", "KRAS", "TP53", "STK11", "KEAP1", "RB1", "ERBB2", "MET")

# ── 0.5 Race/ethnicity colour palette (colour-blind safe) ────────────────────
RACE_COLORS <- c(
  "White (non-Hispanic)"   = "#4477AA",
  "Asian"                  = "#EE6677",
  "Black/African American" = "#228833",
  "Hispanic/Latino"        = "#CCBB44",
  "Other/Unknown"          = "#BBBBBB"
)

RACE_LEVELS <- names(RACE_COLORS)

message("✓ Setup complete. Data directory: ", DATA_DIR)
