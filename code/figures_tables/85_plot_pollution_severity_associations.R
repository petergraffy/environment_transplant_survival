#!/usr/bin/env Rscript

runtime_source <- file.path("code", "r_runtime.R")
if (!file.exists(runtime_source)) runtime_source <- "r_runtime.R"
if (file.exists(runtime_source)) {
  source(runtime_source)
  ensure_user_library()
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(scales)
  library(tibble)
})

baseline_path <- file.path(
  "output",
  "prior_year_pollution_baseline_severity_associations",
  "prior_year_pollution_baseline_severity_associations.csv"
)
tv_path <- file.path(
  "output",
  "timevarying_pollution_severity_models_with_listing_year",
  "timevarying_pollution_severity_results.csv"
)
fig_dir <- file.path("output", "figures", "pollution_severity_associations")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

baseline_results <- read_csv(baseline_path, show_col_types = FALSE)
tv_results <- read_csv(tv_path, show_col_types = FALSE)

continuous_baseline <- baseline_results %>%
  filter(measure == "mean_difference") %>%
  mutate(model_type = "baseline") %>%
  transmute(
    model_type,
    organ_label,
    outcome_label,
    pollutant,
    exposure,
    n,
    estimate,
    conf_low,
    conf_high,
    p_value,
    measure
  )

continuous_tv <- tv_results %>%
  transmute(
    model_type,
    organ_label,
    outcome_label,
    pollutant,
    exposure,
    n = candidate_organ_episodes,
    estimate,
    conf_low,
    conf_high,
    p_value,
    measure
  )

plot_dat <- bind_rows(continuous_baseline, continuous_tv) %>%
  mutate(
    model_type = recode(model_type, baseline = "A. Baseline severity at listing", time_varying = "B. Time-updated severity during follow-up"),
    model_type = factor(model_type, levels = c("A. Baseline severity at listing", "B. Time-updated severity during follow-up")),
    pollutant = factor(pollutant, levels = c("pm25", "o3", "no2")),
    outcome_panel = case_when(
      organ_label == "Kidney" ~ "Kidney: dialysis duration",
      organ_label == "Liver" ~ "Liver: MELD/PELD",
      organ_label == "Heart" ~ "Heart: US-CRS proxy",
      organ_label == "Lung" ~ "Lung: urgency proxy",
      TRUE ~ outcome_label
    ),
    outcome_panel = factor(outcome_panel, levels = rev(c(
      "Kidney: dialysis duration",
      "Liver: MELD/PELD",
      "Heart: US-CRS proxy",
      "Lung: urgency proxy"
    )))
  )

pollutant_labels <- c(
  pm25 = 'PM[2.5]~"per 5 ug/m"^3',
  o3 = 'O[3]~"per 10 ppb"',
  no2 = 'NO[2]~"per 10 ppb"'
)
pollutant_colors <- c(pm25 = "#0072B2", o3 = "#009E73", no2 = "#D55E00")
dodge <- position_dodge(width = 0.72)

base_theme <- theme_minimal(base_size = 18) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(color = "grey15"),
    axis.title.x = element_text(margin = margin(t = 8)),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.justification = "center",
    legend.text = element_text(size = 16),
    legend.title = element_blank(),
    legend.margin = margin(t = 6)
  )

plot_continuous_panel <- function(dat, x_label) {
  ggplot(dat, aes(x = estimate, y = outcome_panel, color = pollutant)) +
  geom_vline(xintercept = 0, color = "grey35", linewidth = 0.45) +
  geom_point(size = 2.8, position = dodge) +
  geom_errorbar(aes(xmin = conf_low, xmax = conf_high), orientation = "y", width = 0.26, linewidth = 0.9, position = dodge) +
  scale_color_manual(values = pollutant_colors, breaks = names(pollutant_labels), labels = parse(text = pollutant_labels)) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
  labs(x = x_label, y = NULL, title = unique(dat$model_type)) +
  base_theme +
  theme(
    plot.title = element_text(face = "bold", size = 20, hjust = 0),
    plot.margin = margin(8, 18, 4, 8)
  )
}

p_baseline <- plot_continuous_panel(
  filter(plot_dat, model_type == "A. Baseline severity at listing"),
  "Adjusted mean difference per prior-year exposure increment"
) +
  theme(legend.position = "none")

p_tv <- plot_continuous_panel(
  filter(plot_dat, model_type == "B. Time-updated severity during follow-up"),
  "Adjusted mean difference per monthly exposure increment"
) +
  theme(axis.text.y = element_blank(), legend.position = "bottom")

combined <- (p_baseline | p_tv) + plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

png_path <- file.path(fig_dir, "figure4_pollution_severity_association_forest.png")
pdf_path <- file.path(fig_dir, "figure4_pollution_severity_association_forest.pdf")
ggsave(png_path, combined, width = 16.5, height = 7.2, dpi = 360, bg = "white")
ggsave(pdf_path, combined, width = 16.5, height = 7.2, device = cairo_pdf, bg = "white")

write_csv(
  tibble(
    figure = "figure4_pollution_severity_association_forest",
    path = normalizePath(c(png_path, pdf_path), winslash = "/", mustWork = FALSE)
  ),
  file.path(fig_dir, "pollution_severity_association_forest_manifest.csv")
)

message("Wrote pollution severity association figure to ", normalizePath(fig_dir, winslash = "/"))
