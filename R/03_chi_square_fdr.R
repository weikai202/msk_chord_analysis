# =============================================================================
# 03_chi_square_fdr.R — Unadjusted chi-square / Fisher tests + BH-FDR
# =============================================================================
# STATISTICAL RATIONALE (统计原理):
#
#   【卡方检验 vs Fisher精确检验】
#     Chi-square test assumes expected cell counts ≥ 5 in all cells of the
#     2×k contingency table (gene mutation status × race). When any expected
#     count < 5 (common for rare mutations in small racial subgroups), we
#     switch to Fisher's exact test, which enumerates the exact hypergeometric
#     distribution without large-sample assumptions.
#     Decision rule implemented: if min(expected) < 5 → Fisher, else chi-square.
#
#   【多重检验校正 — Benjamini-Hochberg FDR】
#     We perform 8 genes × 1 test each = 8 simultaneous tests.
#     The familywise error rate (FWER) controlled by Bonferroni is overly
#     conservative for exploratory genomic analyses. Instead, BH-FDR controls
#     the expected proportion of false discoveries among rejected hypotheses
#     at level q = 0.05. BH procedure:
#       1. Sort p-values: p_(1) ≤ p_(2) ≤ ... ≤ p_(m)
#       2. Find largest k such that p_(k) ≤ (k/m) × q
#       3. Reject H₀ for all tests i ≤ k
#     This is less conservative than Bonferroni while still providing
#     meaningful false discovery control.
#
#   【效应量 — Cramér's V】
#     A significant p-value alone does not indicate practical importance,
#     especially with large N. Cramér's V measures association strength:
#       V = √(χ² / (N × (min(r,c) - 1)))
#     V ≈ 0.1 small, 0.3 medium, 0.5 large (Cohen, 1988).
# =============================================================================

source(here::here("R", "00_setup.R"))
nsclc_df <- readRDS(file.path(DATA_DIR, "nsclc_complete_case.rds"))

# ── 3.1 Test function ─────────────────────────────────────────────────────────
test_gene_race <- function(df, gene) {
  mut_col <- paste0("MUT_", gene)
  # Drop missing race
  sub_df  <- df %>%
    filter(!is.na(RACE_ETH), !is.na(.data[[mut_col]])) %>%
    mutate(MUT = .data[[mut_col]])

  tbl <- table(sub_df$RACE_ETH, sub_df$MUT)

  # Choose test based on expected counts
  expected <- chisq.test(tbl)$expected
  use_fisher <- any(expected < 5)

  if (use_fisher) {
    res <- fisher.test(tbl, simulate.p.value = TRUE, B = 10000)
    test_used <- "Fisher (simulated)"
    stat <- NA_real_
  } else {
    res <- chisq.test(tbl)
    test_used <- "Chi-square"
    stat <- res$statistic
  }

  # Cramér's V
  n    <- sum(tbl)
  k    <- min(dim(tbl)) - 1
  cramerV <- sqrt(chisq.test(tbl, simulate.p.value = use_fisher,
                              B = 10000)$statistic / (n * k))

  # Prevalence per race group
  prev_by_race <- sub_df %>%
    group_by(RACE_ETH) %>%
    summarise(prev = mean(MUT, na.rm = TRUE) * 100, .groups = "drop") %>%
    pivot_wider(names_from = RACE_ETH, values_from = prev,
                names_prefix = "prev_pct_")

  tibble(
    Gene          = gene,
    n_total       = n,
    test_used     = test_used,
    statistic     = as.numeric(stat),
    p_value       = res$p.value,
    cramers_v     = as.numeric(cramerV)
  ) %>%
    bind_cols(prev_by_race)
}

# ── 3.2 Run all gene tests ────────────────────────────────────────────────────
results_raw <- map_dfr(DRIVER_GENES, ~ test_gene_race(nsclc_df, .x))

# ── 3.3 BH-FDR correction ────────────────────────────────────────────────────
results_fdr <- results_raw %>%
  mutate(
    q_value     = p.adjust(p_value, method = "BH"),
    sig_q05     = q_value < 0.05,
    sig_label   = case_when(
      q_value < 0.001 ~ "***",
      q_value < 0.01  ~ "**",
      q_value < 0.05  ~ "*",
      TRUE            ~ "ns"
    )
  ) %>%
  arrange(q_value)

print(results_fdr %>% select(Gene, n_total, test_used, statistic,
                               p_value, q_value, sig_label, cramers_v))

write_csv(results_fdr, file.path(OUTPUT_TAB, "03_chisq_fdr_nsclc.csv"))

# ── 3.4 Bubble plot: effect size vs -log10(q) ─────────────────────────────────
# STATISTICAL RATIONALE (统计原理):
#   Volcano-style plot plots effect magnitude (Cramér's V on x-axis) against
#   statistical significance (-log10(q) on y-axis). Genes in the upper-right
#   quadrant are both statistically significant and clinically meaningful.

p_bubble <- results_fdr %>%
  ggplot(aes(x = cramers_v, y = -log10(q_value + 1e-10),
             size = n_total, colour = sig_q05, label = Gene)) +
  geom_point(alpha = 0.8) +
  geom_text(vjust = -0.8, size = 3.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             colour = "red", linewidth = 0.5) +
  scale_colour_manual(values = c("TRUE" = "#CC3311", "FALSE" = "#AAAAAA"),
                      labels = c("TRUE" = "q < 0.05", "FALSE" = "q ≥ 0.05"),
                      name   = "") +
  scale_size_continuous(name = "N samples", range = c(3, 10)) +
  labs(
    title    = "Gene–race association: effect size vs. significance (BH-FDR)",
    subtitle = "Dashed line: q = 0.05 threshold | NSCLC cohort",
    x        = "Cramér's V (effect size)",
    y        = expression(-log[10](q~value))
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "right")

ggsave(file.path(OUTPUT_FIG, "03_bubble_chisq_fdr.pdf"),
       p_bubble, width = 7, height = 5, dpi = 300)

message("✓ Chi-square + FDR analysis complete.")
