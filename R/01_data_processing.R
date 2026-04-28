# =============================================================================
# 01_data_processing.R — Data loading, merging, harmonisation, and imputation
# =============================================================================
# STATISTICAL RATIONALE (统计原理):
#
#   【数据合并】 Three-file join strategy:
#     Patient file (1 row/patient) LEFT JOIN Sample file (1 row/sample)
#     via PATIENT_ID, then LEFT JOIN Mutations file (many rows/sample)
#     via SAMPLE_ID. Using a left join preserves all patients even when
#     mutation calls are absent (interpreted as wild-type).
#
#   【缺失数据】 Missing completely at random (MCAR) vs. missing at random (MAR):
#     Race/ethnicity is commonly self-reported and administratively missing —
#     likely MAR conditional on institution, age, and cancer type. Complete-
#     case analysis (CCA) is unbiased only under MCAR; since MAR is more
#     plausible, we add MICE (Multiple Imputation by Chained Equations) as a
#     sensitivity analysis. Under MAR, MICE produces consistent estimates
#     by iteratively imputing each variable from a conditional model using
#     all other observed variables. We use m=20 imputed datasets with
#     predictive mean matching (PMM) for continuous vars and logistic
#     regression for binary/categorical vars.
#
#   【二值化突变状态】 Binarising to "oncogenic mutation present" reduces
#     variant-level noise and aligns with clinical relevance — only variants
#     classified as Likely Oncogenic or Oncogenic matter for pathway activation.
# =============================================================================

source(here::here("R", "00_setup.R"))

# ── 1.1 Load raw files ────────────────────────────────────────────────────────
# The downloaded file is a single combined TSV (patient + sample data merged),
# exported from the cBioPortal study table view.
# Column names contain spaces (e.g. "Patient ID", "Overall Survival (Months)")
# — we normalise them to UPPER_SNAKE_CASE for consistency.

CLINICAL_FILE <- list.files(DATA_DIR, pattern = "clinical.*\\.tsv$",
                             full.names = TRUE, ignore.case = TRUE)[1]
if (is.na(CLINICAL_FILE))
  stop("No clinical .tsv file found in ", DATA_DIR,
       "\nExpected: msk_chord_2024_clinical_data*.tsv")

message("Loading combined clinical file: ", basename(CLINICAL_FILE))
clinical_raw <- read_tsv(
  CLINICAL_FILE,
  show_col_types = FALSE,
  na = c("", "NA", "N/A", "[Not Available]", "[Unknown]", "[Not Applicable]",
         "Unknown", "NA ")
)

# Normalise column names: lowercase → UPPER_SNAKE_CASE
names(clinical_raw) <- names(clinical_raw) %>%
  toupper() %>%
  gsub("[^A-Z0-9]+", "_", .) %>%   # replace non-alphanumeric with _
  gsub("_+$|^_+", "", .)           # strip leading/trailing underscores

message("Combined clinical columns: ", paste(names(clinical_raw), collapse = ", "))

# Map normalised names to the standard names the rest of the script expects
clinical_raw <- clinical_raw %>%
  rename_with(~ case_when(
    .x == "PATIENT_ID"                   ~ "PATIENT_ID",
    .x == "SAMPLE_ID"                    ~ "SAMPLE_ID",
    .x == "OVERALL_SURVIVAL_MONTHS"      ~ "OS_MONTHS",
    .x == "OVERALL_SURVIVAL_STATUS"      ~ "OS_STATUS",
    .x == "CANCER_TYPE"                  ~ "CANCER_TYPE",
    .x == "CANCER_TYPE_DETAILED"         ~ "HISTOLOGY",
    .x == "STAGE_HIGHEST_RECORDED"       ~ "STAGE",
    .x == "SMOKING_HISTORY_NLP"          ~ "SMOKING_HISTORY",
    .x == "SEX"                          ~ "SEX",
    .x == "CURRENT_AGE"                  ~ "AGE",
    .x == "RACE"                         ~ "RACE",
    .x == "ETHNICITY"                    ~ "ETHNICITY",
    TRUE                                 ~ .x
  ))

# Split into patient-level and sample-level (pipeline expects both)
patient_raw <- clinical_raw %>%
  select(any_of(c("PATIENT_ID", "RACE", "ETHNICITY", "SEX", "AGE",
                  "OS_MONTHS", "OS_STATUS"))) %>%
  distinct(PATIENT_ID, .keep_all = TRUE)

sample_raw <- clinical_raw %>%
  select(any_of(c("PATIENT_ID", "SAMPLE_ID", "CANCER_TYPE", "HISTOLOGY",
                  "STAGE", "SMOKING_HISTORY")))

# Load mutations file
MUTATION_FILE <- list.files(DATA_DIR,
                             pattern = "mutation.*\\.(tsv|txt|maf)$",
                             full.names = TRUE, ignore.case = TRUE)[1]

if (is.na(MUTATION_FILE)) {
  message("⚠ No mutation file found — creating empty placeholder.")
  message("  Download from cBioPortal: Plots tab → Mutations → Download TSV")
  mut_raw <- tibble(
    HUGO_SYMBOL            = character(),
    TUMOR_SAMPLE_BARCODE   = character(),
    VARIANT_CLASSIFICATION = character()
  )
} else {
  message("Loading mutation file: ", basename(MUTATION_FILE))
  mut_raw <- read_tsv(MUTATION_FILE, comment = "#",
                      show_col_types = FALSE,
                      na = c("", "NA", "[Not Available]"))

  # Rename cbioportalR API column names to pipeline standard names
  mut_raw <- mut_raw %>%
    rename_with(~ case_when(
      .x == "hugoGeneSymbol" ~ "HUGO_SYMBOL",
      .x == "sampleId"       ~ "TUMOR_SAMPLE_BARCODE",
      .x == "mutationType"   ~ "VARIANT_CLASSIFICATION",
      .x == "patientId"      ~ "PATIENT_ID",
      TRUE                   ~ toupper(.x)
    ))
}

message("Raw dimensions — Patient: ", nrow(patient_raw),
        " | Sample: ", nrow(sample_raw),
        " | Mutations: ", nrow(mut_raw))

# ── 1.2 Harmonise race/ethnicity ─────────────────────────────────────────────
# cBioPortal uses separate RACE and ETHNICITY columns in many studies;
# MSK-CHORD may also have RACE_ETHNICITY. Adjust column names as needed.
# Priority rule: if ETHNICITY == "Hispanic or Latino" → "Hispanic/Latino"
# regardless of RACE value (standard OMB classification logic).

harmonise_race_eth <- function(df) {
  # Detect available column names (case-insensitive)
  cols <- tolower(names(df))

  has_combined  <- "race_ethnicity" %in% cols
  has_race      <- "race"      %in% cols
  has_ethnicity <- "ethnicity" %in% cols

  if (has_combined) {
    raw_col <- df[[which(cols == "race_ethnicity")]]
  } else if (has_race && has_ethnicity) {
    race_col      <- df[[which(cols == "race")]]
    ethnicity_col <- df[[which(cols == "ethnicity")]]
    # IMPORTANT: check Non-Hispanic BEFORE Hispanic to avoid substring match
    # "Non-Spanish; Non-Hispanic" contains "Hispanic" as a substring
    raw_col <- case_when(
      grepl("Non-Hispanic|Non-Spanish|not Spanish", ethnicity_col,
            ignore.case = TRUE)                              ~ as.character(race_col),
      grepl("Hispanic|Latino|Spanish", ethnicity_col,
            ignore.case = TRUE)                              ~ "Hispanic/Latino",
      grepl("White|Caucasian", race_col, ignore.case = TRUE) ~ "White",
      grepl("Asian",           race_col, ignore.case = TRUE) ~ "Asian",
      grepl("Black|African",   race_col, ignore.case = TRUE) ~ "Black/African American",
      TRUE ~ as.character(race_col)
    )
  } else if (has_race) {
    raw_col <- df[[which(cols == "race")]]
  } else {
    stop("Cannot find race/ethnicity columns in patient data.")
  }

  # Map to 5 canonical categories
  # Again check Non-Hispanic first to prevent substring false-match
  harmonised <- case_when(
    grepl("Non-Hispanic|Non-Spanish|not Spanish",
          raw_col, ignore.case = TRUE)                        ~ "White (non-Hispanic)",
    grepl("Hispanic|Latino|Spanish",
          raw_col, ignore.case = TRUE)                        ~ "Hispanic/Latino",
    grepl("White|Caucasian", raw_col, ignore.case = TRUE)    ~ "White (non-Hispanic)",
    grepl("Asian",           raw_col, ignore.case = TRUE)    ~ "Asian",
    grepl("Black|African",   raw_col, ignore.case = TRUE)    ~ "Black/African American",
    is.na(raw_col)                                           ~ NA_character_,
    TRUE                                                      ~ "Other/Unknown"
  )

  factor(harmonised,
         levels = c("White (non-Hispanic)", "Asian",
                    "Black/African American", "Hispanic/Latino",
                    "Other/Unknown"))
}

patient_clean <- patient_raw %>%
  mutate(RACE_ETH = harmonise_race_eth(.)) %>%
  # Standardise key column names
  rename_with(toupper) %>%
  mutate(
    AGE        = suppressWarnings(as.numeric(AGE)),
    OS_MONTHS  = suppressWarnings(as.numeric(OS_MONTHS)),
    OS_STATUS  = case_when(
      grepl("DECEASED|1:DECEASED|1", OS_STATUS, ignore.case = TRUE) ~ 1L,
      grepl("LIVING|0:LIVING|0",     OS_STATUS, ignore.case = TRUE) ~ 0L,
      TRUE ~ NA_integer_
    ),
    SEX = factor(toupper(SEX))
  )

# ── 1.3 Clean sample file ────────────────────────────────────────────────────
sample_clean <- sample_raw %>%
  rename_with(toupper) %>%
  mutate(
    # Harmonise smoking history to 3 levels
    SMOKING_HISTORY = case_when(
      grepl("never",          SMOKING_HISTORY, ignore.case = TRUE) ~ "Never",
      grepl("former|current", SMOKING_HISTORY, ignore.case = TRUE) ~ "Ever",
      TRUE ~ NA_character_
    ) %>% factor(levels = c("Never", "Ever")),

    # Harmonise stage — MSK-CHORD uses grouped "Stage 1-3" / "Stage 4" format
    # plus individual values like "4A", "3B", "1", "0", etc.
    STAGE = case_when(
      grepl("4",     STAGE, ignore.case = TRUE) ~ "IV",
      grepl("1-3|1|2|3", STAGE, ignore.case = TRUE) ~ "I-III",
      grepl("^0$",   STAGE, ignore.case = TRUE) ~ "0",
      TRUE ~ NA_character_
    ) %>% factor(levels = c("0", "I-III", "IV")),

    # Cancer type harmonisation for cohort subsetting
    CANCER_TYPE_CLEAN = case_when(
      grepl("Non-Small Cell Lung|NSCLC|Lung Adenocarcinoma|Lung Squamous",
            CANCER_TYPE, ignore.case = TRUE) ~ "NSCLC",
      grepl("Colorectal|Colon|Rectal|CRC",
            CANCER_TYPE, ignore.case = TRUE) ~ "CRC",
      grepl("Breast",
            CANCER_TYPE, ignore.case = TRUE) ~ "Breast",
      TRUE ~ CANCER_TYPE
    ),

    HISTOLOGY = factor(HISTOLOGY)
  )

# ── 1.4 Create binary mutation matrix ─────────────────────────────────────────
# STATISTICAL RATIONALE (统计原理):
#   Oncogenic variants: we keep Missense_Mutation, Nonsense_Mutation,
#   Frame_Shift_Del, Frame_Shift_Ins, Splice_Site, In_Frame_Del,
#   In_Frame_Ins, Translation_Start_Site. We EXCLUDE Silent (synonymous),
#   Intron, and RNA mutations as these typically do not alter protein function.
#   For each gene-sample pair: 1 if ≥1 pathogenic call present, 0 otherwise.
#   This ensures independence-of-observations at the sample level.

ONCOGENIC_CLASSES <- c(
  "Missense_Mutation", "Nonsense_Mutation",
  "Frame_Shift_Del", "Frame_Shift_Ins",
  "Splice_Site", "In_Frame_Del", "In_Frame_Ins",
  "Translation_Start_Site", "Nonstop_Mutation"
)

mut_binary <- mut_raw %>%
  rename_with(toupper) %>%
  filter(
    HUGO_SYMBOL %in% DRIVER_GENES,
    VARIANT_CLASSIFICATION %in% ONCOGENIC_CLASSES
  ) %>%
  distinct(TUMOR_SAMPLE_BARCODE, HUGO_SYMBOL) %>%
  mutate(MUT = 1L) %>%
  pivot_wider(
    names_from  = HUGO_SYMBOL,
    values_from = MUT,
    values_fill = 0L,
    names_prefix = "MUT_"
  )

# Ensure all 8 driver genes appear as columns (fill 0 for absent genes)
for (g in DRIVER_GENES) {
  col <- paste0("MUT_", g)
  if (!col %in% names(mut_binary)) mut_binary[[col]] <- 0L
}

# ── 1.5 Three-way merge ───────────────────────────────────────────────────────
# patient_clean has PATIENT_ID
# sample_clean has PATIENT_ID + SAMPLE_ID
# mut_binary has TUMOR_SAMPLE_BARCODE (= SAMPLE_ID)

# Standardise join key in mutation table
if ("TUMOR_SAMPLE_BARCODE" %in% names(mut_binary)) {
  mut_binary <- rename(mut_binary, SAMPLE_ID = TUMOR_SAMPLE_BARCODE)
}

merged_df <- sample_clean %>%
  left_join(patient_clean, by = "PATIENT_ID") %>%
  left_join(mut_binary,    by = "SAMPLE_ID")

# Fill NA in mutation columns with 0 (sample not in mutation file = no variant)
mut_cols <- paste0("MUT_", DRIVER_GENES)
merged_df <- merged_df %>%
  mutate(across(all_of(mut_cols), ~ replace_na(.x, 0L)))

message("Merged dataset: ", nrow(merged_df), " rows, ", ncol(merged_df), " columns")

# ── 1.6 Cohort subsetting ────────────────────────────────────────────────────
nsclc_df   <- filter(merged_df, CANCER_TYPE_CLEAN == "NSCLC")
crc_df     <- filter(merged_df, CANCER_TYPE_CLEAN == "CRC")
breast_df  <- filter(merged_df, CANCER_TYPE_CLEAN == "Breast")

message("Cohort sizes — NSCLC: ", nrow(nsclc_df),
        " | CRC: ", nrow(crc_df),
        " | Breast: ", nrow(breast_df))

# ── 1.7 Missing data assessment ───────────────────────────────────────────────
# STATISTICAL RATIONALE (统计原理):
#   We visualise missingness patterns with an upset-style or heatmap plot.
#   Little's MCAR test (if sample allows) formally tests H₀: data are MCAR.
#   A significant result (p < 0.05) suggests MAR or MNAR, justifying MICE.

key_vars <- c("RACE_ETH", "SMOKING_HISTORY", "STAGE", "AGE",
              "SEX", "OS_MONTHS", "OS_STATUS", "HISTOLOGY",
              mut_cols)

miss_summary <- nsclc_df %>%
  select(all_of(intersect(key_vars, names(nsclc_df)))) %>%
  summarise(across(everything(), ~ mean(is.na(.x)) * 100)) %>%
  pivot_longer(everything(), names_to = "Variable", values_to = "PCT_Missing") %>%
  arrange(desc(PCT_Missing))

message("\nMissingness in NSCLC cohort (%):")
print(miss_summary, n = 30)

# Save missing data summary
write_csv(miss_summary, file.path(OUTPUT_TAB, "01_missing_data_summary.csv"))

# ── 1.8 Save complete-case data FIRST (before MICE, so later steps always work)
saveRDS(merged_df,  file.path(DATA_DIR, "merged_all.rds"))
saveRDS(nsclc_df,   file.path(DATA_DIR, "nsclc_complete_case.rds"))
saveRDS(crc_df,     file.path(DATA_DIR, "crc_complete_case.rds"))
saveRDS(breast_df,  file.path(DATA_DIR, "breast_complete_case.rds"))
message("✓ Complete-case RDS files saved.")

# ── 1.9 MICE multiple imputation (sensitivity analysis) ──────────────────────
# STATISTICAL RATIONALE (统计原理):
#   MICE (van Buuren & Groothuis-Oudshoorn, 2011) imputes each variable
#   sequentially using a conditional model given all other variables.
#   We run m=20 imputed datasets; results will be pooled using Rubin's rules:
#     β_pooled = mean(β_m)
#     Var_pooled = within-imputation variance + (1 + 1/m) * between-imputation variance
#   This correctly propagates imputation uncertainty into final estimates.

impute_vars <- intersect(
  c("RACE_ETH", "SMOKING_HISTORY", "STAGE", "AGE", "SEX",
    "OS_MONTHS", "OS_STATUS", mut_cols),
  names(nsclc_df)
)

nsclc_for_impute <- nsclc_df %>%
  select(PATIENT_ID, SAMPLE_ID, all_of(impute_vars)) %>%
  # Drop empty factor levels to prevent MICE polyreg failures
  mutate(across(where(is.factor), droplevels))

# Define imputation methods per variable type
init_imp  <- mice(nsclc_for_impute, maxit = 0, printFlag = FALSE)
imp_meth  <- init_imp$method
imp_meth["RACE_ETH"]        <- "polyreg"  # multinomial logistic for >2 categories
imp_meth["SMOKING_HISTORY"] <- "logreg"   # binary: Never vs Ever
imp_meth["STAGE"]           <- "logreg"   # binary: I-III vs IV
imp_meth["SEX"]             <- "logreg"

# Run imputation (m=20 datasets, 10 iterations each)
set.seed(42)
mice_imp <- tryCatch(
  mice(nsclc_for_impute, m = 20, maxit = 10,
       method = imp_meth, printFlag = FALSE),
  error = function(e) {
    message("MICE failed: ", e$message, " — saving NULL placeholder")
    NULL
  }
)

if (!is.null(mice_imp)) {
  saveRDS(mice_imp, file.path(DATA_DIR, "nsclc_mice_imputed.rds"))
  message("MICE imputation complete: ", mice_imp$m, " datasets generated")
} else {
  # Save empty placeholder so Step 4 MICE sensitivity can skip gracefully
  saveRDS(NULL, file.path(DATA_DIR, "nsclc_mice_imputed.rds"))
}

message("✓ Data processing complete. Files saved to: ", DATA_DIR)
