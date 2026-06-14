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
  library(readr)
  library(scales)
  library(stringr)
  library(tibble)
})

in_path <- file.path(
  "output",
  "primary_waitlist_period_pollution_cox",
  "primary_waitlist_period_pollution_cox_results.csv"
)
fig_dir <- file.path("output", "figures", "primary_waitlist_period_pollution_cox")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

theme_primary <- function(base_size = 10.5) {
  theme_minimal(base_size = base_size) +
    theme(
      text = element_text(color = "grey10"),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(linewidth = 0.25, color = "grey86"),
      axis.title.y = element_blank(),
      axis.text = element_text(color = "grey20"),
      plot.title = element_text(face = "bold", size = rel(1.18), margin = margin(b = 3)),
      plot.subtitle = element_text(color = "grey25", margin = margin(b = 8)),
      plot.caption = element_text(color = "grey35", hjust = 0, margin = margin(t = 8)),
      strip.text = element_text(face = "bold", color = "grey10"),
      strip.background = element_rect(fill = "grey94", color = "grey72", linewidth = 0.45),
      legend.position = "top",
      legend.title = element_blank(),
      legend.justification = "left",
      legend.text = element_text(size = rel(0.95)),
      plot.margin = margin(8, 12, 8, 8)
    )
}

pollutant_levels <- c("pm25", "o3", "no2")
pollutant_labels <- c(
  pm25 = "PM₂.₅ per 5 µg/m³",
  o3 = "O₃ per 10 ppb",
  no2 = "NO₂ per 10 ppb"
)
organ_levels <- c("HR", "KI", "LI", "LU")
organ_labels <- c(HR = "Heart", KI = "Kidney", LI = "Liver", LU = "Lung")
pollutant_colors <- c(pm25 = "#0072B2", o3 = "#009E73", no2 = "#D55E00")
organ_colors <- c(Heart = "#0072B2", Kidney = "#D55E00", Liver = "#009E73", Lung = "#CC79A7")
pollutant_shapes <- c(pm25 = 16, o3 = 17, no2 = 15)

format_ci <- function(hr, low, high) {
  sprintf("%.2f (%.2f-%.2f)", hr, low, high)
}

format_p <- function(p) {
  case_when(
    is.na(p) ~ "",
    p < 0.001 ~ "<0.001",
    TRUE ~ sprintf("%.3f", p)
  )
}

results <- read_csv(in_path, show_col_types = FALSE) %>%
  mutate(
    organ = factor(organ, levels = organ_levels, labels = organ_labels[organ_levels]),
    pollutant = factor(exposure, levels = pollutant_levels, labels = pollutant_labels[pollutant_levels]),
    pollutant_key = factor(exposure, levels = pollutant_levels),
    organ_key = factor(as.character(recode(as.character(organ), !!!setNames(organ_levels, organ_labels[organ_levels]))), levels = organ_levels),
    hr_ci = format_ci(hazard_ratio, conf_low, conf_high),
    p_label = format_p(p_value),
    n_label = comma(n),
    event_label = comma(adverse_events),
    time_label = if_else(exposure_data_end_year == 2025, "2006-2025", "2006-2023")
  )

caption_text <- paste(
  "Cause-specific Cox models for death or deterioration delisting.",
  "Models adjust for age, sex, race, listing year, baseline organ score, and stratified baseline hazard by listing center.",
  "PM2.5 and ozone follow-up are censored at 2023-12-31; NO2 follow-up is censored at 2025-12-31."
)
short_caption <- "Adjusted cause-specific Cox models; points are HRs and bars are 95% CIs."

single_axis_dat <- results %>%
  mutate(
    organ = factor(organ, levels = rev(organ_labels[organ_levels])),
    pollutant = factor(pollutant, levels = pollutant_labels[pollutant_levels]),
    pollutant_key = factor(as.character(pollutant_key), levels = pollutant_levels)
  )

p_single_axis <- ggplot(
  single_axis_dat,
  aes(x = hazard_ratio, y = organ, color = organ, shape = pollutant_key)
) +
  geom_vline(xintercept = 1, linewidth = 0.35, color = "grey35") +
  geom_errorbar(
    aes(xmin = conf_low, xmax = conf_high),
    orientation = "y",
    width = 0,
    linewidth = 0.55,
    position = position_dodge(width = 0.58)
  ) +
  geom_point(size = 2.7, stroke = 0.8, position = position_dodge(width = 0.58)) +
  scale_color_manual(values = organ_colors) +
  scale_shape_manual(values = pollutant_shapes, breaks = pollutant_levels, labels = pollutant_labels) +
  scale_x_log10(
    breaks = c(0.8, 1.0, 1.25, 1.5, 2.0),
    labels = label_number(accuracy = 0.01),
    limits = c(0.78, 2.02),
    expand = expansion(mult = c(0.01, 0.03))
  ) +
  guides(
    color = "none",
    shape = guide_legend(order = 1, override.aes = list(color = "grey20", size = 3))
  ) +
  labs(
    title = "Waitlist-Period Air Pollution and Adverse Waitlist Outcomes",
    subtitle = "Cause-specific hazard ratios for death or deterioration delisting",
    x = "Hazard ratio, log scale",
    caption = short_caption
  ) +
  theme_primary(base_size = 10.8) +
  theme(
    legend.box = "vertical",
    legend.spacing.y = unit(2, "pt"),
    legend.margin = margin(0, 0, 2, 0)
  )

forest_by_organ <- results %>%
  mutate(
    pollutant = factor(pollutant, levels = rev(pollutant_labels[pollutant_levels])),
    label_x = conf_high * 1.035
  )

p_by_organ <- ggplot(forest_by_organ, aes(x = hazard_ratio, y = pollutant, color = pollutant_key)) +
  geom_vline(xintercept = 1, linewidth = 0.35, color = "grey35") +
  geom_errorbar(aes(xmin = conf_low, xmax = conf_high), orientation = "y", width = 0.18, linewidth = 0.55) +
  geom_point(size = 2.2) +
  geom_text(
    aes(x = label_x, label = hr_ci),
    color = "grey10",
    hjust = 0,
    size = 3.0,
    show.legend = FALSE
  ) +
  facet_wrap(~organ, ncol = 2) +
  scale_color_manual(values = pollutant_colors, breaks = pollutant_levels, labels = pollutant_labels) +
  scale_x_log10(
    breaks = c(0.8, 1.0, 1.25, 1.5, 2.0),
    labels = label_number(accuracy = 0.01),
    limits = c(0.78, 2.75),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  coord_cartesian(clip = "off") +
  labs(
    title = "Waitlist-Period Air Pollution and Adverse Waitlist Outcomes",
    subtitle = NULL,
    x = "Cause-specific hazard ratio, log scale",
    caption = NULL
  ) +
  theme_primary(base_size = 10.5) +
  theme(
    legend.position = "none",
    panel.spacing = unit(10, "pt"),
    strip.background = element_rect(fill = "grey92", color = "grey65", linewidth = 0.6),
    panel.border = element_rect(fill = NA, color = "grey70", linewidth = 0.55)
  )

forest_by_pollutant <- results %>%
  mutate(
    organ = factor(organ, levels = rev(organ_labels[organ_levels])),
    label_x = 2.22
  )

p_by_pollutant <- ggplot(forest_by_pollutant, aes(x = hazard_ratio, y = organ, color = pollutant_key)) +
  geom_vline(xintercept = 1, linewidth = 0.35, color = "grey35") +
  geom_errorbar(aes(xmin = conf_low, xmax = conf_high), orientation = "y", width = 0.18, linewidth = 0.55) +
  geom_point(size = 2.2) +
  geom_text(
    aes(x = label_x, label = hr_ci),
    color = "grey10",
    hjust = 0,
    size = 3.0,
    show.legend = FALSE
  ) +
  facet_wrap(~pollutant, ncol = 1) +
  scale_color_manual(values = pollutant_colors, breaks = pollutant_levels, labels = pollutant_labels) +
  scale_x_log10(
    breaks = c(0.8, 1.0, 1.25, 1.5, 2.0),
    labels = label_number(accuracy = 0.01),
    limits = c(0.78, 3.0),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  coord_cartesian(clip = "off") +
  labs(
    title = "Primary Cause-Specific Cox Models by Pollutant",
    subtitle = "Hazard ratios for death or deterioration delisting across transplant organs",
    x = "Cause-specific hazard ratio, log scale",
    caption = short_caption
  ) +
  theme_primary(base_size = 10.5) +
  theme(legend.position = "none")

heatmap_dat <- results %>%
  mutate(
    organ = factor(organ, levels = rev(organ_labels[organ_levels])),
    pollutant = factor(pollutant, levels = pollutant_labels[pollutant_levels]),
    heat_label = sprintf("%.2f", hazard_ratio)
  )

p_heatmap <- ggplot(heatmap_dat, aes(x = pollutant, y = organ, fill = hazard_ratio)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = heat_label), size = 3.4, color = "grey10") +
  scale_fill_gradient2(
    low = "#3B7EA1",
    mid = "white",
    high = "#B14A3B",
    midpoint = 1,
    limits = c(0.9, 1.9),
    oob = squish,
    breaks = c(1.0, 1.25, 1.5, 1.75),
    name = "HR"
  ) +
  labs(
    title = "Primary Model Hazard Ratios",
    subtitle = "Cause-specific hazard ratio for adverse waitlist outcome",
    x = NULL,
    y = NULL,
    caption = short_caption
  ) +
  theme_primary(base_size = 10.5) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "right"
  )

save_plot <- function(plot, stem, width, height) {
  png_path <- file.path(fig_dir, paste0(stem, ".png"))
  pdf_path <- file.path(fig_dir, paste0(stem, ".pdf"))
  ggsave(png_path, plot, width = width, height = height, dpi = 600, bg = "white")
  ggsave(pdf_path, plot, width = width, height = height, device = grDevices::pdf, bg = "white")
  c(png = png_path, pdf = pdf_path)
}

paths <- c(
  save_plot(p_single_axis, "primary_waitlist_pollution_cox_forest_single_axis", 8.4, 5.2),
  save_plot(p_by_organ, "primary_waitlist_pollution_cox_forest_by_organ", 8.8, 6.4),
  save_plot(p_by_pollutant, "primary_waitlist_pollution_cox_forest_by_pollutant", 8.8, 7.2),
  save_plot(p_heatmap, "primary_waitlist_pollution_cox_hr_heatmap", 7.4, 4.8)
)

plot_table <- results %>%
  transmute(
    organ = as.character(organ),
    pollutant = as.character(pollutant),
    followup_window = time_label,
    n,
    people,
    centers,
    adverse_events,
    hazard_ratio,
    conf_low,
    conf_high,
    p_value,
    hr_ci,
    p_label
  )

write_csv(plot_table, file.path(fig_dir, "primary_waitlist_pollution_cox_plot_table.csv"))
writeLines(
  c(
    "Figure caption:",
    caption_text,
    "",
    "Primary exposure definition:",
    "Observed waitlist-period average pollution at candidate ZCTA. PM2.5 and ozone are evaluated through 2023; NO2 through 2025.",
    "",
    "Model:",
    "Surv(followup_days_pollutant, adverse_event) ~ waitlist_period_pollution + age + sex + race + listing_year + organ_score + strata(listing_center)."
  ),
  file.path(fig_dir, "primary_waitlist_pollution_cox_caption.txt")
)
write_csv(
  tibble(figure_path = unname(paths)),
  file.path(fig_dir, "primary_waitlist_pollution_cox_figure_manifest.csv")
)

message("Wrote primary waitlist-period pollution Cox figures to ", fig_dir)
