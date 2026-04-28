# =============================================================================
# 08_export_web_figures.R — Convert PDF figures to PNG for GitHub Pages website
# =============================================================================
# Sources:
#   1. PDF figures from Desktop/figures/ (or output/figures/ as fallback)
#      converted to PNG via the magick package.
#   2. Two additional figures (bubble plot, forest OR) generated from
#      saved CSV tables since those PDFs were not produced by the pipeline.
# Run AFTER 07_run_all.R has completed.
# =============================================================================

# ── Locate project root ───────────────────────────────────────────────────────
# here::here() resolves from RStudio's working directory, which may not be the
# project folder.  We search upward from this script's own path first; if that
# fails (e.g. interactive console paste) we fall back to a fixed path.
.locate_root <- function() {
  # 1. Try the path of the currently-sourced file (works when source()'d)
  src <- tryCatch(normalizePath(sys.frame(1)$ofile, mustWork = FALSE),
                  error = function(e) "")
  if (nchar(src) > 0) {
    # script is in <root>/R/ — go up one level
    candidate <- dirname(dirname(src))
    if (file.exists(file.path(candidate, "R", "00_setup.R")))
      return(candidate)
  }
  # 2. Try here::here() in case the working directory is already correct
  h <- tryCatch(here::here(), error = function(e) "")
  if (nchar(h) > 0 && file.exists(file.path(h, "R", "00_setup.R")))
    return(h)
  # 3. Walk upward from getwd()
  parts <- strsplit(normalizePath(getwd(), mustWork = FALSE), "[/\\\\]")[[1]]
  for (i in seq(length(parts), 1)) {
    candidate <- paste(parts[seq_len(i)], collapse = "/")
    if (file.exists(file.path(candidate, "R", "00_setup.R")))
      return(candidate)
  }
  # 4. Hard-coded fallback for this machine
  fallback <- "C:/Users/kai/msk_chord_analysis"
  if (file.exists(file.path(fallback, "R", "00_setup.R"))) return(fallback)
  stop("Cannot find project root. Run  setwd('C:/Users/kai/msk_chord_analysis')  then try again.")
}

PROJ_ROOT <- .locate_root()
message("Project root: ", PROJ_ROOT)
setwd(PROJ_ROOT)           # ensures here::here() works for the rest of the script

library(here)
source(file.path(PROJ_ROOT, "R", "00_setup.R"))

# ── Paths ─────────────────────────────────────────────────────────────────────
# Prefer Desktop/figures if it exists; fall back to output/figures
PDF_DIR <- if (dir.exists(file.path(Sys.getenv("USERPROFILE"), "Desktop", "figures"))) {
  file.path(Sys.getenv("USERPROFILE"), "Desktop", "figures")
} else {
  here::here("output", "figures")
}
message("Reading PDFs from: ", PDF_DIR)

WEB_FIG <- here::here("docs", "figures")
if (!dir.exists(WEB_FIG)) dir.create(WEB_FIG, recursive = TRUE)

# ── Install magick if needed ──────────────────────────────────────────────────
if (!"magick" %in% installed.packages()[, "Package"])
  install.packages("magick", repos = "https://cloud.r-project.org")
library(magick)

# ── Helper: convert one PDF page to PNG ──────────────────────────────────────
pdf_to_png <- function(pdf_name, out_name, density = 180) {
  pdf_path <- file.path(PDF_DIR, pdf_name)
  if (!file.exists(pdf_path)) {
    message("  ⚠ not found, skipping: ", pdf_name)
    return(invisible(NULL))
  }
  img <- image_read_pdf(pdf_path, density = density)
  # KM plots from survminer are multi-page; keep only page 1
  img <- img[1]
  img <- image_background(img, "white")
  image_write(img, file.path(WEB_FIG, out_name), format = "png")
  message("  ✓ ", pdf_name, " → ", out_name)
}

message("\nConverting PDFs to PNG ...")
pdf_to_png("02_mutation_prevalence_all_cohorts.pdf", "fig_01_prevalence_all.png",  density = 180)
pdf_to_png("02_mutation_prevalence_crc.pdf",         "fig_02_prevalence_crc.png",  density = 180)
pdf_to_png("02_mutation_prevalence_breast.pdf",      "fig_03_prevalence_breast.png", density = 180)
pdf_to_png("05_km_race_egfr_nsclc.pdf",              "fig_04_km_race_egfr.png",    density = 180)
pdf_to_png("05_km_egfr_only.pdf",                    "fig_05_km_egfr_only.png",    density = 180)
pdf_to_png("06_replication_forest_all_cohorts.pdf",  "fig_06_replication.png",     density = 180)
pdf_to_png("04_sensitivity_cc_vs_mice.pdf",          "fig_07_sensitivity.png",     density = 180)

# ── Figure 8: Bubble plot (generated from CSV) ────────────────────────────────
message("\nGenerating bubble plot from CSV ...")
chisq_df <- read_csv(
  here::here("output", "tables", "03_chisq_fdr_nsclc.csv"),
  show_col_types = FALSE
)

use_repel <- "ggrepel" %in% installed.packages()[, "Package"]
if (!use_repel) install.packages("ggrepel", repos = "https://cloud.r-project.org")
library(ggrepel)

fig_bubble <- ggplot(chisq_df,
  aes(x = cramers_v, y = -log10(q_value + 1e-10),
      size = n_total, colour = sig_q05, label = Gene)) +
  geom_point(alpha = 0.85) +
  geom_text_repel(size = 4.5, fontface = "bold", show.legend = FALSE,
                   box.padding = 0.4) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             colour = "#CC3311", linewidth = 0.6) +
  annotate("text",
           x = max(chisq_df$cramers_v, na.rm = TRUE) * 0.88,
           y = -log10(0.05) + 2.5,
           label = "q = 0.05", colour = "#CC3311", size = 4) +
  scale_colour_manual(
    values = c("TRUE" = "#CC3311", "FALSE" = "#AAAAAA"),
    labels = c("TRUE" = "q < 0.05", "FALSE" = "q ≥ 0.05"),
    name = "BH-FDR"
  ) +
  scale_size_continuous(name = "N samples", range = c(4, 12)) +
  labs(
    title    = "Gene–race association: effect size vs. significance",
    subtitle = "Cramér's V (x) vs −log₁₀(q) (y) | BH-FDR corrected",
    x = "Cramér's V (effect size)",
    y = expression(-log[10](q~value)),
    caption = "MSK-CHORD 2024 | NSCLC n = 7,775"
  ) +
  theme_bw(base_size = 13) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))

ggsave(file.path(WEB_FIG, "fig_08_bubble.png"),
       fig_bubble, width = 9, height = 6, dpi = 150, bg = "white")
message("  ✓ fig_08_bubble.png")

# ── Figure 9: Adjusted OR forest plot (generated from CSV) ───────────────────
message("Generating forest OR plot from CSV ...")
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
      p_adj < 0.001 ~ "***", p_adj < 0.01 ~ "**",
      p_adj < 0.05 ~ "*",   TRUE ~ ""
    )
  )

fig_forest <- ggplot(or_df,
  aes(y = Race_Group, x = OR, xmin = CI_low, xmax = CI_high,
      colour = Race_Group)) +
  geom_vline(xintercept = 1, linetype = "dashed",
             colour = "grey40", linewidth = 0.6) +
  geom_errorbarh(height = 0.28, linewidth = 0.8) +
  geom_point(size = 4, shape = 18) +
  geom_text(aes(x = pmin(CI_high + 0.08, 6.5), label = sig_label),
            hjust = 0, size = 4, colour = "black") +
  scale_colour_manual(values = RACE_COLORS[-1], guide = "none") +
  scale_x_log10(
    breaks = c(0.25, 0.5, 1, 2, 4),
    labels = c("0.25", "0.5", "1", "2", "4"),
    limits = c(0.18, 8)
  ) +
  facet_wrap(~ Gene, ncol = 2) +
  labs(
    title    = "Adjusted odds ratios for driver gene mutations by race/ethnicity",
    subtitle = "Reference: White (non-Hispanic)  |  Adjusted for: smoking, stage, histology, sex, age",
    x = "Odds ratio (log scale)", y = NULL,
    caption = "* FDR-adjusted p < 0.05  ** p < 0.01  *** p < 0.001  |  MSK-CHORD 2024"
  ) +
  theme_bw(base_size = 13) +
  theme(
    strip.background = element_rect(fill = "#2C3E50"),
    strip.text       = element_text(face = "bold", colour = "white", size = 12),
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold")
  )

ggsave(file.path(WEB_FIG, "fig_09_forest_or.png"),
       fig_forest, width = 10, height = 7, dpi = 150, bg = "white")
message("  ✓ fig_09_forest_or.png")

message("\n✓ All web figures saved to: ", WEB_FIG)
message("  Total PNG files: ",
        length(list.files(WEB_FIG, pattern = "\\.png$")))
