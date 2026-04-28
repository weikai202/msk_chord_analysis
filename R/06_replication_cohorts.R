# =============================================================================
# 06_replication_cohorts.R — Replicate key findings in CRC and Breast cohorts
# =============================================================================
# STATISTICAL RATIONALE (统计原理):
#
#   【复制研究的重要性】
#     Replication in independent cohorts tests whether findings are
#     cancer-type-specific biological effects vs. broader racial disparities
#     in mutation frequencies. Consistent effects across NSCLC, CRC, and
#     Breast suggest a systematic racial/ethnic pattern in tumour genomics
#     (potentially reflecting population-level germline variation or
#     environmental exposure differences), whereas discordant effects
#     suggest histotype-specific biology (e.g., EGFR enrichment in Asian
#     NSCLC is adenocarcinoma-specific, not a pan-cancer phenomenon).
#
#   【多癌种分析注意事项】
#     Each cancer type has different dominant driver genes:
#       NSCLC:  EGFR, KRAS, TP53, STK11, KEAP1
#       CRC:    KRAS, TP53, RB1 (APC not in our gene list)
#       Breast: ERBB2, TP53, RB1
#     Genes irrelevant to a cancer type will show near-zero prevalence —
#     this is expected, not a data error.
# =============================================================================

source(here::here("R", "00_setup.R"))
crc_df    <- readRDS(file.path(DATA_DIR, "crc_complete_case.rds"))
breast_df <- readRDS(file.path(DATA_DIR, "breast_complete_case.rds"))

# ── 6.1 Logistic regression function (reusable) ───────────────────────────────
run_logistic_cohort <- function(df, cohort_name, genes = DRIVER_GENES) {
  df <- df %>%
    filter(!is.na(RACE_ETH)) %>%
    mutate(RACE_ETH = relevel(factor(RACE_ETH), ref = "White (non-Hispanic)"))

  base_cov <- "RACE_ETH + STAGE + SEX + AGE"
  if ("SMOKING_HISTORY" %in% names(df) && mean(!is.na(df$SMOKING_HISTORY)) > 0.5)
    base_cov <- paste0(base_cov, " + SMOKING_HISTORY")

  results <- map_dfr(genes, function(gene) {
    mut_col <- paste0("MUT_", gene)
    if (!mut_col %in% names(df)) return(NULL)

    # Skip genes with <10 events
    if (sum(df[[mut_col]], na.rm = TRUE) < 10) return(NULL)

    formula_str <- paste0(mut_col, " ~ ", base_cov)
    fit <- tryCatch(
      glm(as.formula(formula_str), data = df,
          family = binomial(link = "logit"), na.action = na.omit),
      error = function(e) NULL
    )
    if (is.null(fit)) return(NULL)

    broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
      filter(grepl("^RACE_ETH", term)) %>%
      mutate(
        Gene       = gene,
        Race_Group = sub("RACE_ETH", "", term),
        OR         = estimate,
        CI_low     = conf.low,
        CI_high    = conf.high,
        Cohort     = cohort_name
      ) %>%
      select(Cohort, Gene, Race_Group, OR, CI_low, CI_high, p.value)
  })
  results
}

message("Running CRC cohort logistic regression...")
or_crc <- run_logistic_cohort(crc_df, "CRC")

message("Running Breast cohort logistic regression...")
or_breast <- run_logistic_cohort(breast_df, "Breast")

or_replication <- bind_rows(or_crc, or_breast)
write_csv(or_replication, file.path(OUTPUT_TAB, "06_OR_replication_cohorts.csv"))

# ── 6.2 Three-cohort comparison forest plot ───────────────────────────────────
# Load NSCLC results for comparison
or_nsclc_cc <- read_csv(file.path(OUTPUT_TAB, "04_OR_logistic_nsclc_cc.csv"),
                         show_col_types = FALSE) %>%
  mutate(Cohort = "NSCLC")

or_all_cohorts <- bind_rows(
  or_nsclc_cc %>% select(Cohort, Gene, Race_Group, OR, CI_low, CI_high, p.value),
  or_replication
) %>%
  mutate(
    Cohort     = factor(Cohort, levels = c("NSCLC", "CRC", "Breast")),
    Race_Group = factor(Race_Group, levels = c(
      "Asian", "Black/African American", "Hispanic/Latino", "Other/Unknown"
    ))
  )

# Focus on genes with sufficient events across cohorts
gene_filter <- or_all_cohorts %>%
  group_by(Gene) %>%
  summarise(any_sig = any(p.value < 0.05, na.rm = TRUE)) %>%
  filter(any_sig) %>%
  pull(Gene)

p_multi_cohort <- or_all_cohorts %>%
  filter(Gene %in% gene_filter, !is.na(OR)) %>%
  ggplot(aes(y = Race_Group, x = OR,
             xmin = CI_low, xmax = CI_high,
             colour = Cohort, shape = Cohort)) +
  geom_vline(xintercept = 1, linetype = "dashed",
             colour = "grey40", linewidth = 0.5) +
  geom_errorbarh(height = 0.25, linewidth = 0.6,
                 position = position_dodge(width = 0.6)) +
  geom_point(size = 2.5, position = position_dodge(width = 0.6)) +
  scale_colour_manual(values = c("NSCLC" = "#4477AA",
                                  "CRC"   = "#EE6677",
                                  "Breast"= "#228833"),
                      name = "Cancer cohort") +
  scale_shape_manual(values = c("NSCLC" = 15, "CRC" = 16, "Breast" = 17),
                     name = "Cancer cohort") +
  scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4),
                labels = c("0.25", "0.5", "1", "2", "4")) +
  facet_wrap(~ Gene, ncol = 2) +
  labs(
    title    = "Replication: adjusted ORs across NSCLC, CRC, and Breast cohorts",
    subtitle = "Reference: White (non-Hispanic) | Each cohort adjusted for available confounders",
    x        = "Odds ratio (log scale)", y = NULL
  ) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom",
        strip.background = element_rect(fill = "#E8E8E8"),
        strip.text = element_text(face = "bold"))

ggsave(file.path(OUTPUT_FIG, "06_replication_forest_all_cohorts.pdf"),
       p_multi_cohort, width = 12, height = 16, dpi = 300)

message("✓ Replication analysis complete.")
