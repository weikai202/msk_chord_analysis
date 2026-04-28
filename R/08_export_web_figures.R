# =============================================================================
# 08_export_web_figures.R — Re-create key figures as PNG for the project website
# =============================================================================
# This script reads the saved CSV tables (output/tables/) and produces
# publication-quality PNG figures in docs/figures/.
# It does NOT require the raw .rds data files — only the pipeline outputs.
# Run this script after 07_run_all.R has completed successfully.
# =============================================================================

source(here::here("R", "00_setup.R"))
library(patchwork)

WEB_FIG <- here::here("docs", "figures")
if (!dir.exists(WEB_FIG)) dir.create(WEB_FIG, recursive = TRUE)

# Helper: save PNG at consistent resolution
save_web <- function(p, filename, width = 10, height = 6) {
  ggsave(
    file.path(WEB_FIG, filename), p,
    width = width, height = height, dpi = 150, bg = "white"
  )
  message("  saved: ", filename)
}

# ── Figure 1: Mutation prevalence bar chart ───────────────────────────────────
prev_df <- read_csv(
  here::here("output", "tables", "02_mutation_prevalence_nsclc.csv"),
  show_col_types = FALSE
) %>%
  mutate(
    RACE_ETH = factor(RACE_ETH, levels = names(RACE_COLORS)),
    Gene     = factor(Gene, levels = DRIVER_GENES)
  )

fig1 <- ggplot(prev_df, aes(x = Gene, y = pct, fill = RACE_ETH)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65) +
  geom_errorbar(
    aes(ymin = lo * 100, ymax = hi * 100),
    position = position_dodge(width = 0.75),
    width = 0.25, linewidth = 0.4
  ) +
  scale_fill_manual(values = RACE_COLORS, name = "Race/Ethnicity") +
  scale_y_continuous(
    labels = scales::percent_format(scale = 1),
    expand = expansion(mult = c(0, 0.06))
  ) +
  labs(
    title    = "Driver gene mutation prevalence by race/ethnicity",
    subtitle = "NSCLC cohort (n = 7,775) | Error bars: 95% Wilson CI",
    x = "Driver gene", y = "Mutation prevalence (%)",
    caption = "MSK-CHORD 2024"
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position  = "bottom",
    axis.text.x      = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold")
  ) +
  guides(fill = guide_legend(nrow = 2))

save_web(fig1, "fig_01_prevalence.png", width = 11, height = 6.5)

# ── Figure 2: Bubble plot — effect size vs significance ──────────────────────
chisq_df <- read_csv(
  here::here("output", "tables", "03_chisq_fdr_nsclc.csv"),
  show_col_types = FALSE
)

fig2 <- ggplot(chisq_df,
               aes(x = cramers_v, y = -log10(q_value + 1e-10),
                   size = n_total, colour = sig_q05, label = Gene)) +
  geom_point(alpha = 0.85) +
  ggrepel::geom_text_repel(size = 4.5, fontface = "bold",
                            show.legend = FALSE,
                            box.padding = 0.4) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             colour = "#CC3311", linewidth = 0.6) +
  annotate("text", x = max(chisq_df$cramers_v, na.rm = TRUE) * 0.85,
           y = -log10(0.05) + 1.5,
           label = "q = 0.05", colour = "#CC3311", size = 3.8) +
  scale_colour_manual(
    values = c("TRUE" = "#CC3311", "FALSE" = "#AAAAAA"),
    labels = c("TRUE" = "q < 0.05", "FALSE" = "q ≥ 0.05"),
    name = "BH-FDR"
  ) +
  scale_size_continuous(name = "N samples", range = c(4, 12),
                        guide = guide_legend(override.aes = list(colour = "grey50"))) +
  labs(
    title    = "Gene–race association: effect size vs. statistical significance",
    subtitle = "Cramér's V measures association strength | BH-FDR correction applied",
    x        = "Cramér's V (effect size)",
    y        = expression(-log[10](q~value))
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position  = "right",
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold")
  )

# Suppress if ggrepel not installed
if (!"ggrepel" %in% installed.packages()[, "Package"]) {
  fig2 <- fig2 + geom_text(vjust = -0.8, size = 4, colour = "black")
}

save_web(fig2, "fig_02_bubble.png", width = 9, height = 6)

# ── Figure 3: Adjusted OR forest plot — EGFR, KRAS, STK11 ───────────────────
or_df <- read_csv(
  here::here("output", "tables", "04_OR_logistic_nsclc_cc.csv"),
  show_col_types = FALSE
) %>%
  filter(Gene %in% c("EGFR", "KRAS", "STK11", "KEAP1")) %>%
  mutate(
    Race_Group = factor(Race_Group, levels = c(
      "Asian", "Black/African American", "Hispanic/Latino", "Other/Unknown"
    )),
    Gene = factor(Gene, levels = c("EGFR", "KRAS", "STK11", "KEAP1")),
    sig_label = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      TRUE          ~ ""
    )
  )

fig3 <- ggplot(or_df,
               aes(y = Race_Group, x = OR, xmin = CI_low, xmax = CI_high,
                   colour = Race_Group)) +
  geom_vline(xintercept = 1, linetype = "dashed",
             colour = "grey40", linewidth = 0.6) +
  geom_errorbarh(height = 0.28, linewidth = 0.8) +
  geom_point(size = 3.5, shape = 18) +
  geom_text(aes(x = CI_high + 0.05, label = sig_label),
            hjust = 0, size = 4, colour = "black") +
  scale_colour_manual(values = RACE_COLORS[-1], guide = "none") +
  scale_x_log10(
    breaks = c(0.25, 0.5, 1, 2, 4),
    labels = c("0.25", "0.5", "1", "2", "4"),
    limits = c(0.2, 7)
  ) +
  facet_wrap(~ Gene, ncol = 2) +
  labs(
    title    = "Adjusted odds ratios for driver gene mutations by race/ethnicity",
    subtitle = "Reference: White (non-Hispanic) | Adjusted for: smoking, stage, histology, sex, age",
    x        = "Odds ratio (log scale)", y = NULL,
    caption  = "* FDR-adjusted p < 0.05  ** p < 0.01  *** p < 0.001"
  ) +
  theme_bw(base_size = 13) +
  theme(
    strip.background = element_rect(fill = "#2C3E50"),
    strip.text       = element_text(face = "bold", colour = "white", size = 12),
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold")
  )

save_web(fig3, "fig_03_forest_or.png", width = 10, height = 7)

# ── Figure 4: Median OS dot plot — race × EGFR status ───────────────────────
os_df <- read_csv(
  here::here("output", "tables", "05_km_median_survival.csv"),
  show_col_types = FALSE
) %>%
  filter(RACE_ETH != "Other/Unknown") %>%
  mutate(
    RACE_ETH = factor(RACE_ETH, levels = names(RACE_COLORS)[1:4]),
    EGFR_MUT_LABEL = factor(EGFR_MUT_LABEL,
                             levels = c("EGFR wild-type", "EGFR mutant"))
  )

fig4 <- ggplot(os_df,
               aes(x = RACE_ETH, y = median_os,
                   colour = RACE_ETH, shape = EGFR_MUT_LABEL,
                   group = EGFR_MUT_LABEL)) +
  geom_line(aes(group = EGFR_MUT_LABEL), colour = "grey70",
            linewidth = 0.8, linetype = "dotted") +
  geom_point(size = 5, stroke = 1.2) +
  geom_text(aes(label = round(median_os, 1)),
            vjust = -1, size = 3.5, colour = "black", fontface = "bold") +
  scale_colour_manual(values = unname(RACE_COLORS[1:4]), guide = "none") +
  scale_shape_manual(values = c("EGFR wild-type" = 16, "EGFR mutant" = 17),
                     name = "EGFR status") +
  scale_y_continuous(limits = c(20, 58),
                     breaks  = seq(25, 55, 10)) +
  labs(
    title    = "Median overall survival by race/ethnicity and EGFR mutation status",
    subtitle = "Triangle = EGFR mutant | Circle = EGFR wild-type | Numbers = median OS (months)",
    x = NULL, y = "Median OS (months)",
    caption = "MSK-CHORD 2024 | NSCLC cohort"
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position  = "top",
    panel.grid.minor = element_blank(),
    axis.text.x      = element_text(angle = 20, hjust = 1, face = "bold"),
    plot.title       = element_text(face = "bold")
  )

save_web(fig4, "fig_04_median_os.png", width = 9, height = 6)

# ── Figure 5: Replication — EGFR and KRAS ORs across three cohorts ───────────
rep_df <- read_csv(
  here::here("output", "tables", "06_OR_replication_cohorts.csv"),
  show_col_types = FALSE
)

nsclc_or <- read_csv(
  here::here("output", "tables", "04_OR_logistic_nsclc_cc.csv"),
  show_col_types = FALSE
) %>%
  mutate(Cohort = "NSCLC") %>%
  select(Cohort, Gene, Race_Group, OR, CI_low, CI_high, p.value)

rep_all <- bind_rows(nsclc_or, rep_df) %>%
  filter(
    Gene %in% c("EGFR", "KRAS", "TP53"),
    Race_Group %in% c("Asian", "Black/African American"),
    !is.na(OR), is.finite(OR),
    CI_low > 0, CI_high < 50
  ) %>%
  mutate(
    Cohort     = factor(Cohort, levels = c("NSCLC", "CRC", "Breast")),
    Race_Group = factor(Race_Group,
                        levels = c("Asian", "Black/African American")),
    Gene       = factor(Gene, levels = c("EGFR", "KRAS", "TP53"))
  )

fig5 <- ggplot(rep_all,
               aes(y = Cohort, x = OR, xmin = CI_low, xmax = CI_high,
                   colour = Cohort, shape = Cohort)) +
  geom_vline(xintercept = 1, linetype = "dashed",
             colour = "grey40", linewidth = 0.5) +
  geom_errorbarh(height = 0.3, linewidth = 0.7,
                 position = position_dodge(width = 0.5)) +
  geom_point(size = 3.5, position = position_dodge(width = 0.5)) +
  scale_colour_manual(
    values = c("NSCLC" = "#4477AA", "CRC" = "#EE6677", "Breast" = "#228833"),
    name = "Cancer cohort"
  ) +
  scale_shape_manual(
    values = c("NSCLC" = 15, "CRC" = 16, "Breast" = 17),
    name = "Cancer cohort"
  ) +
  scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4),
                labels  = c("0.25", "0.5", "1", "2", "4")) +
  facet_grid(Race_Group ~ Gene) +
  labs(
    title    = "Replication: adjusted ORs across NSCLC, CRC, and Breast cohorts",
    subtitle = "Reference: White (non-Hispanic) | Consistent EGFR enrichment in Asian NSCLC only",
    x        = "Odds ratio (log scale)", y = NULL,
    caption  = "MSK-CHORD 2024"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position  = "bottom",
    strip.background = element_rect(fill = "#2C3E50"),
    strip.text       = element_text(face = "bold", colour = "white"),
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold")
  )

save_web(fig5, "fig_05_replication.png", width = 11, height = 7)

message("\n✓ All web figures exported to: ", WEB_FIG)
