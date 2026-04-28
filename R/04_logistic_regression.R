# =============================================================================
# 04_logistic_regression.R — Multivariable logistic regression + forest plots
# =============================================================================
# STATISTICAL RATIONALE (统计原理):
#
#   【为什么用逻辑回归】
#     The outcome is binary (mutation present/absent), so ordinary linear
#     regression is inappropriate — predictions can exceed [0,1] and residuals
#     are non-normal. Logistic regression models the log-odds:
#       logit(P(Y=1)) = β₀ + β₁·RACE_ETH + β₂·smoking + β₃·stage + ...
#     The exponentiated coefficient exp(β₁) is the odds ratio (OR): the ratio
#     of odds of mutation for a given racial group vs. the reference
#     (White non-Hispanic), holding all other covariates constant.
#
#   【混杂因素控制 — DAG推导】
#     Based on the directed acyclic graph (DAG):
#       RACE → [smoking, SES, access] → cancer_type → mutation_freq
#       RACE → mutation_freq (direct biological path, e.g. EGFR in Asian NSCLC)
#     Covariates: smoking_history, stage, histologic_subtype, sex, age block
#     back-door paths through these variables.
#     NOTE: We do NOT condition on variables on the causal pathway
#     (e.g., treatment) to avoid collider bias.
#
#   【OR解释】
#     OR > 1: higher odds of mutation in that group vs. reference.
#     OR = 1: no association after adjustment.
#     95% CI not crossing 1 (i.e., p < 0.05) → statistically significant.
#     Clinically meaningful threshold: OR > 1.5 or < 0.67 (arbitrary but
#     commonly used in genomic epidemiology).
#
#   【模型假设检验】
#     Logistic regression requires: no perfect multicollinearity (check VIF),
#     no complete separation (Firth penalised regression as fallback),
#     and sufficient events per variable (EPV ≥ 10 recommended;
#     EPV = min(events, non-events) / number of predictors).
# =============================================================================

source(here::here("R", "00_setup.R"))
nsclc_df  <- readRDS(file.path(DATA_DIR, "nsclc_complete_case.rds"))

mice_path <- file.path(DATA_DIR, "nsclc_mice_imputed.rds")
mice_imp  <- if (file.exists(mice_path)) readRDS(mice_path) else NULL
if (is.null(mice_imp)) message("⚠ MICE file not found — skipping sensitivity analysis")

# ── 4.1 Covariate formula builder ────────────────────────────────────────────
# Collapse HISTOLOGY to 3 categories to avoid quasi-complete separation
# (hundreds of rare subtypes each predict mutation with probability 0 or 1)
nsclc_df <- nsclc_df %>%
  mutate(
    RACE_ETH        = relevel(RACE_ETH, ref = "White (non-Hispanic)"),
    SMOKING_HISTORY = relevel(SMOKING_HISTORY, ref = "Never"),
    STAGE           = relevel(droplevels(STAGE), ref = "I-III"),
    HIST3 = factor(case_when(
      grepl("Adenocarcinoma|LUAD", HISTOLOGY, ignore.case = TRUE) ~ "Adenocarcinoma",
      grepl("Squamous|LUSC",       HISTOLOGY, ignore.case = TRUE) ~ "Squamous",
      TRUE                                                         ~ "Other/NOS"
    ), levels = c("Adenocarcinoma", "Squamous", "Other/NOS"))
  )

base_covars <- "RACE_ETH + SMOKING_HISTORY + STAGE + SEX + AGE + HIST3"

# ── 4.2 Fit one logistic model per gene ──────────────────────────────────────
fit_logistic <- function(df, gene) {
  outcome <- paste0("MUT_", gene)
  if (!outcome %in% names(df)) return(NULL)

  # Check events per variable (EPV)
  events   <- sum(df[[outcome]], na.rm = TRUE)
  n_params <- length(strsplit(base_covars, "\\+")[[1]]) + 4  # rough expansion
  epv      <- events / n_params
  message(sprintf("  Gene: %-8s | events: %d | EPV ≈ %.1f", gene, events, epv))

  formula_str <- paste0(outcome, " ~ ", base_covars)
  fit <- tryCatch(
    glm(as.formula(formula_str), data = df, family = binomial(link = "logit"),
        na.action = na.omit),
    error = function(e) {
      message("  Standard GLM failed for ", gene, ": ", e$message,
              " — trying Firth penalised regression")
      if ("logistf" %in% installed.packages()[,"Package"]) {
        library(logistf)
        logistf(as.formula(formula_str), data = df)
      } else NULL
    }
  )
  fit
}

message("\nFitting logistic models (NSCLC, complete-case):")
models_cc <- setNames(
  map(DRIVER_GENES, ~ fit_logistic(nsclc_df, .x)),
  DRIVER_GENES
)

# ── 4.3 Extract OR + 95% CI with broom::tidy ─────────────────────────────────
extract_or <- function(fit, gene) {
  if (is.null(fit)) return(NULL)

  broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(grepl("^RACE_ETH", term)) %>%
    mutate(
      Gene       = gene,
      Race_Group = sub("RACE_ETH", "", term),
      OR         = estimate,
      CI_low     = conf.low,
      CI_high    = conf.high
    ) %>%
    select(Gene, Race_Group, OR, CI_low, CI_high, p.value)
}

or_table_cc <- map2_dfr(models_cc, DRIVER_GENES,
                         ~ extract_or(.x, .y)) %>%
  mutate(
    p_adj       = p.adjust(p.value, method = "BH"),
    sig_label   = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      TRUE          ~ ""
    ),
    Analysis = "Complete-case"
  )

write_csv(or_table_cc, file.path(OUTPUT_TAB, "04_OR_logistic_nsclc_cc.csv"))

# ── 4.4 MICE sensitivity analysis — pool results with Rubin's rules ───────────
# STATISTICAL RATIONALE (统计原理):
#   After MICE, we fit the same logistic model on each of the m=20 imputed
#   datasets, then pool coefficients using with() + pool() from the mice package.
#   Rubin's rules combine:
#     - Within-imputation variance: average of the m estimated variances
#     - Between-imputation variance: variance of the m point estimates
#   The pooled statistic follows a t-distribution with Barnard-Rubin df.

fit_pooled_logistic <- function(gene) {
  outcome     <- paste0("MUT_", gene)
  formula_str <- paste0(outcome, " ~ ", base_covars)
  fit_all <- with(mice_imp,
    glm(as.formula(formula_str), family = binomial(link = "logit")))
  pooled  <- pool(fit_all)
  summary(pooled, conf.int = TRUE, exponentiate = TRUE) %>%
    as_tibble() %>%
    filter(grepl("^RACE_ETH", term)) %>%
    mutate(
      Gene       = gene,
      Race_Group = sub("RACE_ETH", "", term),
      OR         = estimate,
      CI_low     = `2.5 %`,
      CI_high    = `97.5 %`,
      Analysis   = "MICE"
    ) %>%
    select(Gene, Race_Group, OR, CI_low, CI_high, p.value, Analysis)
}

if (!is.null(mice_imp)) {
  message("\nFitting pooled logistic models (MICE imputed datasets):")
  or_table_mice <- map_dfr(DRIVER_GENES, fit_pooled_logistic)
  write_csv(or_table_mice, file.path(OUTPUT_TAB, "04_OR_logistic_nsclc_mice.csv"))
} else {
  message("Skipping MICE sensitivity analysis (no imputed data).")
  or_table_mice <- tibble(Gene=character(), Race_Group=character(),
                           OR=numeric(), CI_low=numeric(), CI_high=numeric(),
                           p.value=numeric(), Analysis=character())
}

# ── 4.5 Combine CC + MICE for sensitivity comparison ──────────────────────────
or_combined <- bind_rows(
  or_table_cc %>% select(Gene, Race_Group, OR, CI_low, CI_high, p.value, Analysis),
  or_table_mice
) %>%
  mutate(
    Race_Group = factor(Race_Group, levels = c(
      "Asian", "Black/African American", "Hispanic/Latino", "Other/Unknown"
    )),
    sig_label = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      TRUE ~ ""
    )
  )

write_csv(or_combined, file.path(OUTPUT_TAB, "04_OR_logistic_combined.csv"))

# ── 4.6 Forest plot — faceted by gene ────────────────────────────────────────
# STATISTICAL RATIONALE (统计原理):
#   Forest plots display point estimates (OR) with horizontal lines for 95% CI.
#   The vertical reference line at OR=1 represents null (no association).
#   Squares are sized proportional to inverse-variance weight (1/SE²) to
#   visually emphasise more precise estimates.

p_forest <- or_combined %>%
  filter(Analysis == "Complete-case") %>%
  ggplot(aes(y = Race_Group, x = OR,
             xmin = CI_low, xmax = CI_high,
             colour = Race_Group, label = sig_label)) +
  geom_vline(xintercept = 1, linetype = "dashed",
             colour = "grey40", linewidth = 0.5) +
  geom_errorbarh(height = 0.25, linewidth = 0.7) +
  geom_point(size = 3, shape = 15) +
  geom_text(aes(x = CI_high + 0.05), hjust = 0, size = 3,
            colour = "black", nudge_x = 0.1) +
  scale_colour_manual(values = RACE_COLORS[-1], guide = "none") +  # exclude ref
  scale_x_log10(
    breaks = c(0.25, 0.5, 1, 2, 4, 8),
    labels = c("0.25", "0.5", "1", "2", "4", "8")
  ) +
  facet_wrap(~ Gene, ncol = 2, scales = "free_x") +
  labs(
    title    = "Adjusted odds ratios for driver gene mutations by race/ethnicity",
    subtitle = "Reference: White (non-Hispanic) | Adjusted for: smoking, stage, histology, sex, age",
    x        = "Odds ratio (log scale)",
    y        = NULL,
    caption  = "MSK-CHORD 2024 | Complete-case NSCLC cohort | * p<0.05, ** p<0.01, *** p<0.001"
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "#E8E8E8"),
    strip.text       = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(OUTPUT_FIG, "04_forest_logistic_nsclc.pdf"),
       p_forest, width = 11, height = 14, dpi = 300)

# ── 4.7 Sensitivity comparison: CC vs MICE side-by-side ──────────────────────
p_sensitivity <- or_combined %>%
  ggplot(aes(y = Race_Group, x = OR,
             xmin = CI_low, xmax = CI_high,
             colour = Analysis, shape = Analysis)) +
  geom_vline(xintercept = 1, linetype = "dashed",
             colour = "grey40", linewidth = 0.5) +
  geom_errorbarh(height = 0.25, linewidth = 0.6,
                 position = position_dodge(width = 0.4)) +
  geom_point(size = 2.5,
             position = position_dodge(width = 0.4)) +
  scale_x_log10() +
  scale_colour_manual(values = c("Complete-case" = "#4477AA",
                                  "MICE" = "#EE6677"),
                      name = "Analysis") +
  facet_wrap(~ Gene, ncol = 2, scales = "free_x") +
  labs(
    title    = "Complete-case vs. MICE sensitivity analysis",
    subtitle = "NSCLC | Adjusted ORs for RACE_ETH",
    x        = "Odds ratio (log scale)", y = NULL
  ) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom",
        strip.background = element_rect(fill = "#E8E8E8"))

ggsave(file.path(OUTPUT_FIG, "04_sensitivity_cc_vs_mice.pdf"),
       p_sensitivity, width = 11, height = 14, dpi = 300)

message("✓ Logistic regression analysis complete.")
