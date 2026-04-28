# MSK-CHORD 2024 — Racial/Ethnic Disparities in Driver Gene Mutations

## Data download

1. Go to https://www.cbioportal.org/study/summary?id=msk_chord_2024
2. Click **Download** → **Compressed study**
3. Extract the `.tar.gz` to `data/raw/`
4. Confirm these files exist:
   - `data/raw/data_clinical_patient.txt`
   - `data/raw/data_clinical_sample.txt`
   - `data/raw/data_mutations.txt`

## Run the pipeline

```r
# Option A: Run from RStudio (recommended)
setwd("path/to/msk_chord_analysis")
source("R/07_run_all.R")

# Option B: From terminal
Rscript R/07_run_all.R
```

## Output files

| File | Description |
|------|-------------|
| `output/tables/01_missing_data_summary.csv` | % missing per variable |
| `output/tables/02_table1_nsclc.csv` | Table 1 stratified by race/ethnicity |
| `output/tables/02_mutation_prevalence_nsclc.csv` | Raw prevalence per gene/race |
| `output/tables/03_chisq_fdr_nsclc.csv` | Chi-square + BH-FDR results |
| `output/tables/04_OR_logistic_nsclc_cc.csv` | Adjusted ORs (complete-case) |
| `output/tables/04_OR_logistic_nsclc_mice.csv` | Adjusted ORs (MICE imputed) |
| `output/tables/05_cox_hr_nsclc.csv` | Cox PH hazard ratios |
| `output/tables/05_cox_ph_test.csv` | Schoenfeld PH test results |
| `output/tables/06_OR_replication_cohorts.csv` | ORs in CRC + Breast cohorts |
| `output/figures/02_mutation_prevalence_nsclc.pdf` | Grouped bar plot NSCLC |
| `output/figures/03_bubble_chisq_fdr.pdf` | Effect size vs. significance bubble |
| `output/figures/04_forest_logistic_nsclc.pdf` | Logistic forest plot (8 genes) |
| `output/figures/04_sensitivity_cc_vs_mice.pdf` | CC vs MICE sensitivity |
| `output/figures/05_km_race_egfr_nsclc.pdf` | KM curves Race × EGFR |
| `output/figures/05_forest_cox_interaction.pdf` | Cox interaction forest plot |
| `output/figures/05_schoenfeld_residuals.pdf` | PH assumption diagnostic |
| `output/figures/06_replication_forest_all_cohorts.pdf` | 3-cohort replication |

## R session info

R >= 4.3.0 required. Key packages: tidyverse, broom, survival, survminer, mice,
gtsummary, patchwork, scales, viridis.
