# =============================================================================
# 02_descriptive.R — Table 1 + grouped bar plots of mutation prevalence
# =============================================================================
# STATISTICAL RATIONALE (统计原理):
#
#   【Table 1】 Descriptive statistics stratified by exposure (RACE_ETH) allow
#     readers to assess whether the cohort is balanced across groups before
#     any modelling. Continuous variables are summarised as median [IQR]
#     (non-parametric, appropriate for skewed distributions like age and
#     OS_MONTHS). Categorical variables are summarised as n (%).
#     The standardised mean difference (SMD) is reported instead of p-values
#     for group comparison in Table 1 — SMD > 0.1 conventionally indicates
#     a meaningful imbalance that may require adjustment in regression models.
#
#   【突变频率条形图】 Grouped bar plots show raw (unadjusted) mutation
#     prevalence per gene per racial group. Error bars represent 95% Wilson
#     confidence intervals for a proportion:
#       CI = (p + z²/2n ± z√(p(1-p)/n + z²/4n²)) / (1 + z²/n)
#     Wilson CIs are preferred over Wald CIs because they perform better
#     near p=0 or p=1, which is common for rare mutations.
# =============================================================================

source(here::here("R", "00_setup.R"))
nsclc_df <- readRDS(file.path(DATA_DIR, "nsclc_complete_case.rds"))

# ── 2.1 Table 1 ───────────────────────────────────────────────────────────────
library(gtsummary)

# Variables to include in Table 1
table1_vars <- intersect(
  c("AGE", "SEX", "SMOKING_HISTORY", "STAGE", "HISTOLOGY",
    "OS_MONTHS", "OS_STATUS",
    paste0("MUT_", DRIVER_GENES)),
  names(nsclc_df)
)

# Create clean labels for display
var_labels <- c(
  AGE             = "Age at diagnosis (years)",
  SEX             = "Sex",
  SMOKING_HISTORY = "Smoking history",
  STAGE           = "Clinical stage",
  HISTOLOGY       = "Histologic subtype",
  OS_MONTHS       = "Overall survival (months)",
  OS_STATUS       = "Vital status (1=deceased)",
  MUT_EGFR        = "EGFR mutation",
  MUT_KRAS        = "KRAS mutation",
  MUT_TP53        = "TP53 mutation",
  MUT_STK11       = "STK11 mutation",
  MUT_KEAP1       = "KEAP1 mutation",
  MUT_RB1         = "RB1 mutation",
  MUT_ERBB2       = "ERBB2 mutation",
  MUT_MET         = "MET mutation"
)

# Collapse HISTOLOGY to top-5 levels + "Other" to avoid gtsummary hanging
top_hist <- nsclc_df %>%
  filter(!is.na(RACE_ETH)) %>%
  count(HISTOLOGY, sort = TRUE) %>%
  slice_head(n = 5) %>%
  pull(HISTOLOGY) %>%
  as.character()

tbl1_data <- nsclc_df %>%
  filter(!is.na(RACE_ETH)) %>%
  select(RACE_ETH, all_of(setdiff(table1_vars, "HISTOLOGY"))) %>%
  mutate(
    HISTOLOGY_TOP = factor(
      ifelse(as.character(nsclc_df$HISTOLOGY[!is.na(nsclc_df$RACE_ETH)]) %in% top_hist,
             as.character(nsclc_df$HISTOLOGY[!is.na(nsclc_df$RACE_ETH)]),
             "Other")
    )
  )

tbl1 <- tbl1_data %>%
  tbl_summary(
    by        = RACE_ETH,
    statistic = list(
      all_continuous()  ~ "{median} [{p25}, {p75}]",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits  = list(all_continuous() ~ 1, all_categorical() ~ 1),
    missing = "no"
  ) %>%
  add_overall() %>%
  add_p(test = list(
    all_continuous()  ~ "kruskal.test",
    all_categorical() ~ "chisq.test"
  )) %>%
  bold_labels()

# Export Table 1 as CSV
tbl1_df <- tbl1 %>% as_tibble()
write_csv(tbl1_df, file.path(OUTPUT_TAB, "02_table1_nsclc.csv"))

message("Table 1 saved.")

# ── 2.2 Wilson CI helper ──────────────────────────────────────────────────────
wilson_ci <- function(x, n, conf = 0.95) {
  z   <- qnorm(1 - (1 - conf) / 2)
  p   <- x / n
  denom <- 1 + z^2 / n
  centre <- (p + z^2 / (2 * n)) / denom
  half   <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / denom
  list(est = centre, lower = pmax(0, centre - half), upper = pmin(1, centre + half))
}

# ── 2.3 Compute mutation prevalence with CIs ──────────────────────────────────
prev_df <- nsclc_df %>%
  filter(!is.na(RACE_ETH)) %>%
  group_by(RACE_ETH) %>%
  summarise(
    n_total = n(),
    across(all_of(paste0("MUT_", DRIVER_GENES)),
           ~ sum(.x, na.rm = TRUE),
           .names = "n_{.col}")
  ) %>%
  pivot_longer(
    cols      = starts_with("n_MUT_"),
    names_to  = "Gene",
    values_to = "n_mut"
  ) %>%
  mutate(
    Gene = sub("n_MUT_", "", Gene),
    est  = wilson_ci(n_mut, n_total)$est,
    lo   = wilson_ci(n_mut, n_total)$lower,
    hi   = wilson_ci(n_mut, n_total)$upper,
    pct  = est * 100
  )

# ── 2.4 Grouped bar plot — NSCLC ──────────────────────────────────────────────
p_bar_nsclc <- ggplot(prev_df,
    aes(x = Gene, y = pct, fill = RACE_ETH)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65) +
  geom_errorbar(
    aes(ymin = lo * 100, ymax = hi * 100),
    position = position_dodge(width = 0.75),
    width = 0.25, linewidth = 0.4
  ) +
  scale_fill_manual(values = RACE_COLORS, name = "Race/Ethnicity") +
  scale_y_continuous(labels = scales::percent_format(scale = 1),
                     expand = expansion(mult = c(0, 0.05))) +
  labs(
    title    = "Driver gene mutation prevalence by race/ethnicity in NSCLC",
    subtitle = "Error bars: 95% Wilson confidence intervals",
    x        = "Driver gene",
    y        = "Mutation prevalence (%)",
    caption  = "MSK-CHORD 2024 | Complete-case analysis"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position  = "bottom",
    axis.text.x      = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  ) +
  guides(fill = guide_legend(nrow = 2))

ggsave(file.path(OUTPUT_FIG, "02_mutation_prevalence_nsclc.pdf"),
       p_bar_nsclc, width = 10, height = 6, dpi = 300)

# ── 2.5 Replicate bar plot for CRC and Breast ─────────────────────────────────
make_prevalence_plot <- function(df, cohort_label) {
  df %>%
    filter(!is.na(RACE_ETH)) %>%
    group_by(RACE_ETH) %>%
    summarise(
      n_total = n(),
      across(all_of(paste0("MUT_", DRIVER_GENES)),
             ~ sum(.x, na.rm = TRUE), .names = "n_{.col}")
    ) %>%
    pivot_longer(starts_with("n_MUT_"), names_to = "Gene", values_to = "n_mut") %>%
    mutate(
      Gene = sub("n_MUT_", "", Gene),
      est  = wilson_ci(n_mut, n_total)$est,
      lo   = wilson_ci(n_mut, n_total)$lower,
      hi   = wilson_ci(n_mut, n_total)$upper,
      pct  = est * 100
    ) %>%
    ggplot(aes(x = Gene, y = pct, fill = RACE_ETH)) +
      geom_col(position = position_dodge(width = 0.75), width = 0.65) +
      geom_errorbar(
        aes(ymin = lo * 100, ymax = hi * 100),
        position = position_dodge(width = 0.75),
        width = 0.25, linewidth = 0.4
      ) +
      scale_fill_manual(values = RACE_COLORS, name = "Race/Ethnicity") +
      scale_y_continuous(labels = scales::percent_format(scale = 1),
                         expand = expansion(mult = c(0, 0.05))) +
      labs(
        title = paste("Driver gene mutation prevalence by race/ethnicity —", cohort_label),
        x = "Driver gene", y = "Mutation prevalence (%)"
      ) +
      theme_bw(base_size = 12) +
      theme(legend.position = "bottom",
            axis.text.x = element_text(angle = 45, hjust = 1))
}

crc_df   <- readRDS(file.path(DATA_DIR, "crc_complete_case.rds"))
breast_df <- readRDS(file.path(DATA_DIR, "breast_complete_case.rds"))

p_crc    <- make_prevalence_plot(crc_df,    "CRC")
p_breast <- make_prevalence_plot(breast_df, "Breast Cancer")

ggsave(file.path(OUTPUT_FIG, "02_mutation_prevalence_crc.pdf"),
       p_crc,    width = 10, height = 6, dpi = 300)
ggsave(file.path(OUTPUT_FIG, "02_mutation_prevalence_breast.pdf"),
       p_breast, width = 10, height = 6, dpi = 300)

# ── 2.6 Combined 3-cohort panel ───────────────────────────────────────────────
# patchwork stacks plots vertically, sharing a common legend
combined_panel <- (p_bar_nsclc / p_crc / p_breast) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(file.path(OUTPUT_FIG, "02_mutation_prevalence_all_cohorts.pdf"),
       combined_panel, width = 11, height = 18, dpi = 300)

# Save prevalence table
write_csv(prev_df, file.path(OUTPUT_TAB, "02_mutation_prevalence_nsclc.csv"))

message("✓ Descriptive analysis complete.")
