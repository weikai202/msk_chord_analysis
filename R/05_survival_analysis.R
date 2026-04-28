# =============================================================================
# 05_survival_analysis.R — KM curves, log-rank test, Cox PH with interaction
# =============================================================================
# STATISTICAL RATIONALE (统计原理):
#
#   【生存分析基础 — 删失数据】
#     Survival analysis handles right-censored observations: patients who
#     have not experienced the event (death) by last follow-up contribute
#     partial information. Standard regression ignores censoring and biases
#     estimates downward. The Kaplan-Meier estimator:
#       Ŝ(t) = ∏_{tᵢ ≤ t} [1 - dᵢ/nᵢ]
#     where dᵢ = deaths at time tᵢ, nᵢ = at-risk patients.
#     KM is non-parametric: no distributional assumption on event times.
#
#   【对数秩检验 (Log-rank test)】
#     H₀: survival functions are identical across groups.
#     The test compares observed vs. expected deaths under H₀ at each event
#     time using a chi-square statistic:
#       χ² = Σ_k (O_k - E_k)² / V_k
#     Log-rank is most powerful for proportional hazards (constant HR over
#     time). We apply it here for its interpretability and familiarity.
#     For ≥3 groups we use the overall test; pairwise tests use Bonferroni.
#
#   【Cox比例风险模型 (Cox PH)】
#     h(t|X) = h₀(t) · exp(β₁X₁ + β₂X₂ + ...)
#     The baseline hazard h₀(t) is left unspecified (semi-parametric).
#     exp(β) is the hazard ratio (HR): the multiplicative change in the
#     instantaneous risk of death for a 1-unit increase in the covariate.
#     The key assumption is PROPORTIONAL HAZARDS: HRs are constant over time.
#     Violation means the effect of a covariate changes with time → add
#     time-varying covariate or stratify.
#
#   【交互项 RACE_ETH * EGFR_mut】
#     The interaction term tests whether the effect of EGFR mutation on
#     survival differs across racial groups (effect modification). A
#     significant interaction (p < 0.05) means the HR for EGFR_mut is not
#     the same in all racial groups — racial disparities in OS partly mediated
#     through differential EGFR mutation biology. Interaction on the
#     multiplicative scale (log HR) is the natural scale for Cox models.
#
#   【Schoenfeld残差检验比例风险假设】
#     Schoenfeld residuals are the difference between observed covariate
#     values and their expected values under PH. If PH holds, residuals
#     are uncorrelated with time. cox.zph() tests:
#       H₀: residuals are independent of time (PH holds)
#     A significant result suggests time-varying effect; consider adding
#     an interaction with log(time) or stratifying on that covariate.
# =============================================================================

source(here::here("R", "00_setup.R"))
nsclc_df <- readRDS(file.path(DATA_DIR, "nsclc_complete_case.rds"))

# ── 5.1 Prepare survival data ─────────────────────────────────────────────────
surv_df <- nsclc_df %>%
  filter(
    !is.na(RACE_ETH),
    !is.na(OS_MONTHS),
    !is.na(OS_STATUS),
    OS_MONTHS >= 0
  ) %>%
  mutate(
    RACE_ETH = relevel(factor(RACE_ETH), ref = "White (non-Hispanic)"),
    STAGE    = relevel(droplevels(STAGE), ref = "I-III"),
    # Create 4-level stratification variable for KM curves
    KM_GROUP = interaction(RACE_ETH, factor(MUT_EGFR), sep = " | EGFR="),
    EGFR_MUT_LABEL = factor(MUT_EGFR, levels = c(0, 1),
                             labels = c("EGFR wild-type", "EGFR mutant"))
  )

message("Survival analysis N: ", nrow(surv_df))

# ── 5.2 Kaplan-Meier — RACE_ETH × EGFR mutation ──────────────────────────────
surv_obj <- Surv(time = surv_df$OS_MONTHS, event = surv_df$OS_STATUS)

km_fit <- survfit(surv_obj ~ RACE_ETH + EGFR_MUT_LABEL, data = surv_df)

# Log-rank test
lr_test <- survdiff(surv_obj ~ RACE_ETH + EGFR_MUT_LABEL, data = surv_df)
lr_pval <- 1 - pchisq(lr_test$chisq, df = length(lr_test$n) - 1)
message(sprintf("Log-rank test (all groups): χ²=%.2f, df=%d, p=%.4f",
                lr_test$chisq, length(lr_test$n) - 1, lr_pval))

# ── 5.3 KM plot with survminer ────────────────────────────────────────────────
# Custom colour palette: 2 shades per racial group (solid=mutant, dashed=WT)
race_base_colours <- unname(RACE_COLORS[1:5])
km_colours <- c(rbind(race_base_colours,
                       scales::alpha(race_base_colours, 0.4)))

km_plot <- ggsurvplot(
  km_fit,
  data           = surv_df,
  palette        = km_colours,
  pval           = TRUE,
  pval.method    = TRUE,
  conf.int       = TRUE,
  conf.int.alpha = 0.1,
  risk.table     = TRUE,
  risk.table.height = 0.28,
  xlab           = "Time (months)",
  ylab           = "Overall survival probability",
  title          = "KM curves: Race/Ethnicity × EGFR mutation status (NSCLC)",
  legend.title   = "Group",
  ggtheme        = theme_bw(base_size = 11),
  risk.table.fontsize = 3,
  surv.median.line   = "hv"
)

# Save KM plot
km_pdf_path <- file.path(OUTPUT_FIG, "05_km_race_egfr_nsclc.pdf")
pdf(km_pdf_path, width = 12, height = 9)
print(km_plot)
dev.off()
message("KM plot saved: ", km_pdf_path)

# ── 5.4 KM by EGFR status only (simpler overview) ────────────────────────────
km_egfr <- survfit(surv_obj ~ EGFR_MUT_LABEL, data = surv_df)
p_km_simple <- ggsurvplot(
  km_egfr,
  data        = surv_df,
  palette     = c("#4477AA", "#EE6677"),
  pval        = TRUE,
  conf.int    = TRUE,
  risk.table  = TRUE,
  xlab        = "Time (months)",
  ylab        = "Overall survival probability",
  title       = "KM curves: EGFR mutation status (NSCLC)",
  ggtheme     = theme_bw(base_size = 11)
)

pdf(file.path(OUTPUT_FIG, "05_km_egfr_only.pdf"), width = 10, height = 7)
print(p_km_simple)
dev.off()

# ── 5.5 Multivariable Cox PH with interaction ─────────────────────────────────
cox_formula <- as.formula(paste0(
  "Surv(OS_MONTHS, OS_STATUS) ~ RACE_ETH * MUT_EGFR + ",
  "SMOKING_HISTORY + STAGE + SEX + AGE",
  ifelse("HISTOLOGY" %in% names(surv_df) &&
           mean(!is.na(surv_df$HISTOLOGY)) > 0.7,
         " + HISTOLOGY", "")
))

cox_fit <- coxph(cox_formula, data = surv_df, ties = "efron", x = TRUE)
# ties = "efron" handles tied event times more accurately than "breslow"

message("\nCox model summary:")
print(summary(cox_fit))

# ── 5.6 Extract HR + 95% CI ───────────────────────────────────────────────────
cox_tidy <- broom::tidy(cox_fit, exponentiate = TRUE, conf.int = TRUE) %>%
  rename(HR = estimate, CI_low = conf.low, CI_high = conf.high) %>%
  mutate(
    sig_label = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      TRUE            ~ ""
    )
  )

write_csv(cox_tidy, file.path(OUTPUT_TAB, "05_cox_hr_nsclc.csv"))
print(cox_tidy %>% select(term, HR, CI_low, CI_high, p.value, sig_label), n = 30)

# ── 5.7 Test proportional hazards assumption (Schoenfeld residuals) ───────────
# STATISTICAL RATIONALE (统计原理):
#   cox.zph() fits a weighted linear regression of Schoenfeld residuals on
#   time. If the slope ≠ 0, the effect changes over time → PH violated.
#   We report the global test and per-covariate tests.

ph_test <- cox.zph(cox_fit, transform = "km")
ph_df   <- tibble(
  covariate  = rownames(ph_test$table),
  chi_sq     = ph_test$table[, "chisq"],
  df         = ph_test$table[, "df"],
  p_value    = ph_test$table[, "p"],
  ph_violated = ph_test$table[, "p"] < 0.05
)

write_csv(ph_df, file.path(OUTPUT_TAB, "05_cox_ph_test.csv"))
message("\nProportional hazards test:")
print(ph_df)

# Plot Schoenfeld residuals for variables that may violate PH
pdf(file.path(OUTPUT_FIG, "05_schoenfeld_residuals.pdf"),
    width = 10, height = 8)
par(mfrow = c(2, 3))
for (i in seq_len(min(6, nrow(ph_test$table) - 1))) {
  plot(ph_test[i],
       main  = paste("Schoenfeld residuals:", rownames(ph_test$table)[i]),
       xlab  = "Time", ylab = "Beta(t)")
  abline(h = coef(cox_fit)[i], lty = 2, col = "red")
}
dev.off()

# ── 5.8 Forest plot for Cox results ──────────────────────────────────────────
# Show only main effects and interaction terms of interest
terms_of_interest <- c(
  grep("^RACE_ETH",   cox_tidy$term, value = TRUE),
  grep("^MUT_EGFR",   cox_tidy$term, value = TRUE),
  grep("RACE_ETH.*MUT_EGFR|MUT_EGFR.*RACE_ETH", cox_tidy$term, value = TRUE)
)

p_cox_forest <- cox_tidy %>%
  filter(term %in% terms_of_interest) %>%
  mutate(
    term_clean = term %>%
      sub("RACE_ETH", "Race: ", .) %>%
      sub("MUT_EGFR", "EGFR mut", .) %>%
      sub(":", " × ", .),
    is_interaction = grepl("×|:", term)
  ) %>%
  ggplot(aes(y = reorder(term_clean, HR),
             x = HR, xmin = CI_low, xmax = CI_high,
             colour = is_interaction, label = sig_label)) +
  geom_vline(xintercept = 1, linetype = "dashed",
             colour = "grey40", linewidth = 0.5) +
  geom_errorbarh(height = 0.3, linewidth = 0.7) +
  geom_point(size = 3.5, shape = 18) +
  geom_text(aes(x = CI_high + 0.02), hjust = 0, size = 3.5,
            colour = "black") +
  scale_colour_manual(
    values = c("FALSE" = "#4477AA", "TRUE" = "#CC3311"),
    labels = c("FALSE" = "Main effect", "TRUE" = "Interaction term"),
    name   = "Term type"
  ) +
  scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4),
                labels = c("0.25", "0.5", "1", "2", "4")) +
  labs(
    title    = "Multivariable Cox PH: HR for RACE_ETH × EGFR interaction",
    subtitle = "Adjusted for smoking, stage, histology, sex, age | Reference: White (non-Hispanic)",
    x        = "Hazard ratio (log scale)", y = NULL,
    caption  = "NSCLC cohort | * p<0.05, ** p<0.01, *** p<0.001"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

ggsave(file.path(OUTPUT_FIG, "05_forest_cox_interaction.pdf"),
       p_cox_forest, width = 10, height = 7, dpi = 300)

# ── 5.9 Median survival table ────────────────────────────────────────────────
km_summary <- surv_df %>%
  filter(!is.na(RACE_ETH)) %>%
  group_by(RACE_ETH, EGFR_MUT_LABEL) %>%
  summarise(
    n          = n(),
    n_events   = sum(OS_STATUS, na.rm = TRUE),
    median_os  = summary(survfit(Surv(OS_MONTHS, OS_STATUS) ~ 1))$table["median"],
    .groups    = "drop"
  )

write_csv(km_summary, file.path(OUTPUT_TAB, "05_km_median_survival.csv"))
print(km_summary)

message("✓ Survival analysis complete.")
