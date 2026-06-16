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
  "primary_waitlist_period_pollution_subgroup_cox",
  "primary_waitlist_period_pollution_subgroup_cox_results.csv"
)
interaction_path <- file.path(
  "output",
  "primary_waitlist_period_pollution_subgroup_cox",
  "primary_waitlist_period_pollution_subgroup_interactions.csv"
)
fig_dir <- file.path("output", "figures", "primary_waitlist_period_pollution_subgroup_cox")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

pollutant_levels <- c("pm25", "o3", "no2")
pollutant_dodge_levels <- c("no2", "o3", "pm25")
pollutant_labels <- c(
  pm25 = "PM\u2082.\u2085 per 5 \u00b5g/m\u00b3",
  o3 = "O\u2083 per 10 ppb",
  no2 = "NO\u2082 per 10 ppb"
)
pollutant_labels_plotmath <- c(
  pm25 = 'PM[2.5]~"per 5 ug/m"^3',
  o3 = 'O[3]~"per 10 ppb"',
  no2 = 'NO[2]~"per 10 ppb"'
)
pollutant_colors <- c(pm25 = "#0072B2", o3 = "#009E73", no2 = "#D55E00")
organ_levels <- c("Heart", "Kidney", "Liver", "Lung")
subgroup_domain_levels <- c("Age group", "Sex", "Race")
dodge_width <- 0.9

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

theme_subgroup_forest <- function(base_size = 16) {
  theme_minimal(base_size = base_size) +
    theme(
      text = element_text(color = "black"),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(linewidth = 0.18, color = "grey91"),
      panel.grid.major.x = element_line(linewidth = 0.25, color = "grey86"),
      axis.title.y = element_blank(),
      axis.text = element_text(color = "black", size = rel(0.95)),
      axis.text.y = element_text(size = rel(0.95)),
      axis.title.x = element_text(size = rel(1.05)),
      plot.title = element_text(face = "bold", size = rel(1.2), margin = margin(b = 8)),
      strip.text = element_text(face = "bold", color = "black", size = rel(1.0), margin = margin(5, 5, 5, 5)),
      strip.background = element_rect(fill = "grey94", color = "grey68", linewidth = 0.55),
      legend.position = "top",
      legend.title = element_blank(),
      legend.justification = "left",
      legend.margin = margin(0, 0, 4, 0),
      panel.spacing.x = unit(10, "pt"),
      panel.spacing.y = unit(10, "pt"),
      panel.border = element_rect(fill = NA, color = "grey75", linewidth = 0.45),
      plot.margin = margin(8, 12, 8, 8)
    )
}

results <- read_csv(in_path, show_col_types = FALSE) %>%
  filter(is.na(error), !(subgroup_label == "Race" & subgroup_level == "Other/Unknown")) %>%
  mutate(
    organ_label = factor(organ_label, levels = organ_levels),
    subgroup_label = factor(subgroup_label, levels = subgroup_domain_levels),
    exposure = factor(exposure, levels = pollutant_dodge_levels),
    pollutant_label = factor(pollutant_labels[as.character(exposure)], levels = pollutant_labels),
    subgroup_level = case_when(
      subgroup_label == "Age group" ~ factor(subgroup_level, levels = rev(c("<18", "18-39", "40-59", "60+"))),
      subgroup_label == "Sex" ~ factor(
        recode(subgroup_level, F = "Female", M = "Male"),
        levels = rev(c("Female", "Male"))
      ),
      TRUE ~ factor(subgroup_level, levels = rev(c("White", "Black", "Asian")))
    ),
    hr_ci = format_ci(hazard_ratio, conf_low, conf_high),
    label_left = organ_label == "Kidney" &
      subgroup_label == "Race" &
      subgroup_level == "Black" &
      as.character(exposure) == "pm25",
    label_x = if_else(label_left, conf_low * 0.97, conf_high * 1.03),
    label_hjust = if_else(label_left, 1, 0)
  )

interaction_labels <- read_csv(interaction_path, show_col_types = FALSE) %>%
  mutate(
    organ_label = factor(organ_label, levels = organ_levels),
    subgroup_label = factor(subgroup_label, levels = subgroup_domain_levels),
    exposure = factor(exposure, levels = pollutant_levels),
    y = Inf,
    label = paste0(pollutant_labels[as.character(exposure)], " p int. ", format_p(interaction_p_value))
  ) %>%
  group_by(organ_label, subgroup_label) %>%
  summarise(
    label = paste(label, collapse = "\n"),
    .groups = "drop"
  ) %>%
  mutate(
    x = 0.58,
    y = Inf
  )

p_subgroup <- ggplot(
  results,
  aes(x = hazard_ratio, y = subgroup_level, color = exposure)
) +
  geom_vline(xintercept = 1, linewidth = 0.35, color = "grey34") +
  geom_errorbar(
    aes(xmin = conf_low, xmax = conf_high),
    orientation = "y",
    width = 0,
    linewidth = 0.45,
    position = position_dodge(width = dodge_width)
  ) +
  geom_point(size = 2.6, stroke = 0.8, shape = 16, position = position_dodge(width = dodge_width)) +
  geom_text(
    aes(x = label_x, label = hr_ci, hjust = label_hjust, group = exposure),
    color = "black",
    size = 3.45,
    position = position_dodge(width = dodge_width),
    show.legend = FALSE
  ) +
  facet_grid(subgroup_label ~ organ_label, scales = "free_y", space = "free_y") +
  scale_color_manual(values = pollutant_colors, breaks = pollutant_levels, labels = parse(text = pollutant_labels_plotmath[pollutant_levels])) +
  scale_x_log10(
    breaks = c(0.35, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0),
    labels = label_number(accuracy = 0.1),
    limits = c(0.35, 6.45),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  coord_cartesian(clip = "off") +
  guides(
    color = guide_legend(order = 1, override.aes = list(size = 3))
  ) +
  labs(
    title = "Subgroup Associations Between Waitlist-Period Air Pollution and Adverse Waitlist Outcomes",
    x = "Cause-specific hazard ratio, log scale"
  ) +
  theme_subgroup_forest(base_size = 16) +
  theme(
    legend.key.width = unit(24, "pt"),
    legend.text = element_text(size = 15),
    axis.title.x = element_text(margin = margin(t = 10))
  )

p_subgroup_annotated <- p_subgroup +
  geom_text(
    data = interaction_labels,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1.1,
    color = "grey28",
    size = 2.0,
    lineheight = 0.95
  )

png_path <- file.path(fig_dir, "primary_waitlist_pollution_subgroup_forest.png")
pdf_path <- file.path(fig_dir, "primary_waitlist_pollution_subgroup_forest.pdf")
annotated_png_path <- file.path(fig_dir, "primary_waitlist_pollution_subgroup_forest_with_interaction_p.png")
annotated_pdf_path <- file.path(fig_dir, "primary_waitlist_pollution_subgroup_forest_with_interaction_p.pdf")
ggsave(png_path, p_subgroup, width = 20, height = 16, dpi = 320, bg = "white")
ggsave(pdf_path, p_subgroup, width = 20, height = 16, device = cairo_pdf, bg = "white")
ggsave(annotated_png_path, p_subgroup_annotated, width = 20, height = 16, dpi = 320, bg = "white")
ggsave(annotated_pdf_path, p_subgroup_annotated, width = 20, height = 16, device = cairo_pdf, bg = "white")

manifest <- tibble(
  file = c(basename(png_path), basename(pdf_path), basename(annotated_png_path), basename(annotated_pdf_path)),
  path = c(
    normalizePath(png_path, winslash = "/", mustWork = FALSE),
    normalizePath(pdf_path, winslash = "/", mustWork = FALSE),
    normalizePath(annotated_png_path, winslash = "/", mustWork = FALSE),
    normalizePath(annotated_pdf_path, winslash = "/", mustWork = FALSE)
  ),
  width_in = 20,
  height_in = 16,
  description = c(
    "Clean subgroup forest plot for primary waitlist-period pollution cause-specific Cox models.",
    "Clean subgroup forest plot for primary waitlist-period pollution cause-specific Cox models.",
    "Subgroup forest plot with exposure-by-subgroup interaction p-values.",
    "Subgroup forest plot with exposure-by-subgroup interaction p-values."
  )
)
write_csv(manifest, file.path(fig_dir, "primary_waitlist_pollution_subgroup_forest_manifest.csv"))

message("Wrote subgroup forest plot to ", fig_dir)
