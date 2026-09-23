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

in_path <- file.path(
  "output",
  "prior_year_and_timevarying_pollution_subgroup_cox",
  "prior_year_and_timevarying_pollution_subgroup_cox_results.csv"
)
fig_dir <- file.path("output", "figures", "prior_year_and_timevarying_pollution_subgroup_cox")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

pollutant_levels <- c("pm25", "no2")
pollutant_dodge_levels <- c("no2", "pm25")
pollutant_labels_plotmath <- c(
  pm25 = 'PM[2.5]~"per 5 ug/m"^3',
  no2 = 'NO[2]~"per 10 ppb"'
)
pollutant_colors <- c(pm25 = "#0072B2", no2 = "#D55E00")
organ_levels <- c("Heart", "Lung")
subgroup_domain_levels <- c("Age group", "Sex", "Race")
dodge_width <- 0.62

format_ci <- function(hr, low, high) {
  sprintf("%.2f (%.2f-%.2f)", hr, low, high)
}

theme_focused_subgroup <- function(base_size = 26) {
  theme_minimal(base_size = base_size) +
    theme(
      text = element_text(color = "black"),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(linewidth = 0.18, color = "grey91"),
      panel.grid.major.x = element_line(linewidth = 0.25, color = "grey86"),
      axis.title.y = element_blank(),
      axis.text = element_text(color = "black"),
      axis.text.y = element_text(size = rel(1.0)),
      axis.title.x = element_text(size = rel(1.04), margin = margin(t = 11)),
      plot.title = element_text(face = "bold", size = rel(1.08), margin = margin(b = 9)),
      strip.text = element_text(face = "bold", color = "black", size = rel(1.02), margin = margin(7, 7, 7, 7)),
      strip.background = element_rect(fill = "white", color = "grey25", linewidth = 0.75),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.justification = "center",
      legend.text = element_text(size = rel(0.98)),
      panel.spacing.x = unit(14, "pt"),
      panel.spacing.y = unit(12, "pt"),
      panel.border = element_rect(fill = NA, color = "grey55", linewidth = 0.45),
      plot.margin = margin(10, 72, 10, 10)
    )
}

plot_one_model <- function(dat, model_value, title_text) {
  plot_dat <- dat %>%
    filter(
      model == model_value,
      is.na(error),
      organ_label %in% organ_levels,
      pollutant %in% pollutant_levels,
      !(subgroup_label == "Race" & subgroup_level == "Other/Unknown")
    ) %>%
    mutate(
      organ_label = factor(organ_label, levels = organ_levels),
      subgroup_label = factor(subgroup_label, levels = subgroup_domain_levels),
      pollutant = factor(pollutant, levels = pollutant_dodge_levels),
      subgroup_level = case_when(
        subgroup_label == "Age group" ~ factor(subgroup_level, levels = rev(c("<18", "18-39", "40-59", "60+"))),
        subgroup_label == "Sex" ~ factor(subgroup_level, levels = rev(c("Female", "Male"))),
        TRUE ~ factor(subgroup_level, levels = rev(c("White", "Black", "Asian")))
      ),
      hr_ci = format_ci(hazard_ratio, conf_low, conf_high),
      label_x = conf_high * 1.01
    )

  ggplot(plot_dat, aes(x = hazard_ratio, y = subgroup_level, color = pollutant)) +
    geom_vline(xintercept = 1, linewidth = 0.45, color = "grey30") +
    geom_errorbar(
      aes(xmin = conf_low, xmax = conf_high),
      orientation = "y",
      width = 0.10,
      linewidth = 0.7,
      position = position_dodge(width = dodge_width)
    ) +
    geom_point(size = 3.8, shape = 16, position = position_dodge(width = dodge_width)) +
    geom_text(
      aes(x = label_x, label = hr_ci, group = pollutant),
      hjust = 0,
      color = "black",
      size = 5.5,
      position = position_dodge(width = dodge_width),
      show.legend = FALSE
    ) +
    facet_grid(subgroup_label ~ organ_label, scales = "free_y", space = "free_y") +
    scale_color_manual(values = pollutant_colors, breaks = pollutant_levels, labels = parse(text = pollutant_labels_plotmath[pollutant_levels])) +
    scale_x_log10(
      breaks = c(0.75, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0),
      labels = label_number(accuracy = 0.01),
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    coord_cartesian(xlim = c(0.62, 3.2), clip = "off") +
    guides(color = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(size = 4.4))) +
    labs(title = title_text, x = "Cause-specific hazard ratio, log scale") +
    theme_focused_subgroup()
}

results <- read_csv(in_path, show_col_types = FALSE)

p_baseline <- plot_one_model(results, "baseline_prior_year", "A. 1-year exposure before listing")
p_tv <- plot_one_model(results, "time_varying", "B. Time-varying waitlist-period exposure")

combined <- p_baseline / p_tv + plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

png_path <- file.path(fig_dir, "heart_lung_pm25_no2_subgroup_forest.png")
pdf_path <- file.path(fig_dir, "heart_lung_pm25_no2_subgroup_forest.pdf")
ggsave(png_path, combined, width = 21, height = 23, dpi = 360, bg = "white")
ggsave(pdf_path, combined, width = 21, height = 23, device = cairo_pdf, bg = "white")

write_csv(
  tibble(
    figure = "heart_lung_pm25_no2_subgroup_forest",
    path = normalizePath(c(png_path, pdf_path), winslash = "/", mustWork = FALSE)
  ),
  file.path(fig_dir, "heart_lung_pm25_no2_subgroup_forest_manifest.csv")
)

message("Wrote focused heart/lung PM2.5/NO2 subgroup forest plot to ", normalizePath(fig_dir, winslash = "/"))
