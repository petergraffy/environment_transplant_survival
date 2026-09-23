#!/usr/bin/env Rscript

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", cmd_args[grepl("^--file=", cmd_args)])
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(file_arg[[1]], winslash = "/", mustWork = FALSE))
} else {
  file.path(getwd(), "code", "figures_tables")
}
repo_root <- normalizePath(file.path(script_dir, "..", ".."), winslash = "/", mustWork = FALSE)
if (dir.exists(repo_root)) setwd(repo_root)

runtime_source <- file.path(repo_root, "code", "r_runtime.R")
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
  "prior_year_initial_listing_cox_model_set",
  "prior_year_initial_listing_cox_model_set_results.csv"
)
fig_dir <- file.path("output", "figures", "prior_year_initial_listing_cox_model_set")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

organ_levels <- c("Heart", "Kidney", "Liver", "Lung")
pollutant_levels <- c("PM2.5", "O3", "NO2")
model_tier_levels <- c(
  "Unadjusted",
  "Adjusted + SDOH",
  "Adjusted + SDOH + organ score at listing"
)
model_tier_labels <- c(
  "Unadjusted",
  "Adjusted + SDOH",
  "Adjusted + SDOH + organ score"
)
pollutant_labels <- c(
  "PM2.5" = "PM[2.5]",
  "O3" = "O[3]",
  "NO2" = "NO[2]"
)
model_tier_colors <- c(
  "Unadjusted" = "#666666",
  "Adjusted + SDOH" = "#0072B2",
  "Adjusted + SDOH + organ score at listing" = "#D55E00"
)

format_ci <- function(hr, low, high) {
  sprintf("%.2f (%.2f-%.2f)", hr, low, high)
}

pollutant_from_term <- function(x) {
  case_when(
    x == "pm25_prior_5ug" ~ "PM2.5",
    x == "o3_prior_10ppb" ~ "O3",
    x == "no2_prior_10ppb" ~ "NO2",
    TRUE ~ x
  )
}

theme_model_tier_forest <- function(base_size = 18) {
  theme_minimal(base_size = base_size) +
    theme(
      text = element_text(color = "grey10"),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(linewidth = 0.25, color = "grey86"),
      axis.title.y = element_blank(),
      axis.text = element_text(color = "grey15"),
      axis.text.y = element_text(size = rel(0.95)),
      axis.title.x = element_text(size = rel(1.0), margin = margin(t = 10)),
      strip.text = element_text(face = "bold", color = "grey10", size = rel(1.0), margin = margin(7, 7, 7, 7)),
      strip.background = element_rect(fill = "white", color = "grey45", linewidth = 0.7),
      panel.border = element_rect(fill = NA, color = "grey70", linewidth = 0.5),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = rel(0.92)),
      legend.key.width = unit(22, "pt"),
      panel.spacing.x = unit(10, "pt"),
      panel.spacing.y = unit(9, "pt"),
      plot.margin = margin(8, 48, 8, 8)
    )
}

prepare_plot_dat <- function(results, model_type_filter = NULL) {
  dat <- results %>%
    mutate(
      organ_label = factor(organ_label, levels = organ_levels),
      pollutant = factor(pollutant_from_term(exposure), levels = pollutant_levels),
      model_tier = factor(model_tier, levels = model_tier_levels, labels = model_tier_labels),
      model_tier_raw = factor(model_tier, levels = model_tier_labels),
      hr_ci = format_ci(hazard_ratio, conf_low, conf_high),
      label_x = conf_high * 1.015
    )

  if (!is.null(model_type_filter)) {
    dat <- dat %>% filter(str_detect(model_type, model_type_filter))
  }

  dat
}

make_single_pollutant_plot <- function(plot_dat) {
  single_plot_dat <- plot_dat %>%
    mutate(
      model_tier = factor(as.character(model_tier), levels = rev(model_tier_labels)),
      label_x = conf_high * 1.012
    )

  ggplot(single_plot_dat, aes(x = hazard_ratio, y = model_tier, color = model_tier)) +
    geom_vline(xintercept = 1, linewidth = 0.45, color = "grey35") +
    geom_errorbar(
      aes(xmin = conf_low, xmax = conf_high),
      orientation = "y",
      width = 0.14,
      linewidth = 0.65,
      show.legend = FALSE
    ) +
    geom_point(size = 3.0, show.legend = FALSE) +
    geom_text(
      aes(x = label_x, label = hr_ci),
      hjust = 0,
      size = 3.7,
      color = "grey10",
      show.legend = FALSE
    ) +
    facet_grid(
      organ_label ~ pollutant,
      labeller = labeller(pollutant = as_labeller(pollutant_labels, label_parsed))
    ) +
    scale_color_manual(values = unname(model_tier_colors[model_tier_levels]), breaks = model_tier_labels) +
    scale_x_log10(
      breaks = c(0.8, 1.0, 1.25, 1.5),
      labels = label_number(accuracy = 0.01),
      expand = expansion(mult = c(0.01, 0.08))
    ) +
    coord_cartesian(xlim = c(0.78, max(1.62, max(single_plot_dat$conf_high, na.rm = TRUE) * 1.24)), clip = "off") +
    labs(x = "Cause-specific hazard ratio, log scale") +
    theme_model_tier_forest(base_size = 16) +
    theme(
      legend.position = "none",
      axis.text.y = element_text(size = rel(0.88)),
      plot.margin = margin(8, 58, 8, 8)
    )
}

make_all_spec_plot <- function(plot_dat) {
  all_spec_dat <- plot_dat %>%
    mutate(
      model_type_display = recode(
        model_type,
        "Single-pollutant PM2.5" = "Single pollutant",
        "Single-pollutant O3" = "Single pollutant",
        "Single-pollutant NO2" = "Single pollutant",
        "Multipollutant PM2.5 + NO2" = "PM2.5 + NO2",
        "Multipollutant PM2.5 + NO2 + O3" = "PM2.5 + NO2 + O3"
      ),
      model_type_display = factor(
        model_type_display,
        levels = rev(c("Single pollutant", "PM2.5 + NO2", "PM2.5 + NO2 + O3"))
      )
    )

  ggplot(all_spec_dat, aes(x = hazard_ratio, y = model_type_display, color = model_tier)) +
    geom_vline(xintercept = 1, linewidth = 0.38, color = "grey35") +
    geom_errorbar(
      aes(xmin = conf_low, xmax = conf_high),
      orientation = "y",
      width = 0.12,
      linewidth = 0.55,
      position = position_dodge(width = 0.72)
    ) +
    geom_point(size = 2.35, position = position_dodge(width = 0.72)) +
    facet_grid(
      organ_label ~ pollutant,
      labeller = labeller(pollutant = as_labeller(pollutant_labels, label_parsed))
    ) +
    scale_color_manual(values = unname(model_tier_colors[model_tier_levels]), breaks = model_tier_labels) +
    scale_x_log10(
      breaks = c(0.8, 1.0, 1.25, 1.5),
      labels = label_number(accuracy = 0.01),
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    coord_cartesian(xlim = c(0.78, max(1.58, max(all_spec_dat$conf_high, na.rm = TRUE) * 1.05)), clip = "off") +
    guides(color = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(size = 3.2))) +
    labs(x = "Cause-specific hazard ratio, log scale") +
    theme_model_tier_forest(base_size = 15) +
    theme(plot.margin = margin(8, 12, 8, 8))
}

results <- read_csv(in_path, show_col_types = FALSE)

single_dat <- prepare_plot_dat(results, "^Single-pollutant")
all_dat <- prepare_plot_dat(results)

single_plot <- make_single_pollutant_plot(single_dat)
all_spec_plot <- make_all_spec_plot(all_dat)

single_png <- file.path(fig_dir, "initial_listing_single_pollutant_model_tier_forest.png")
single_pdf <- file.path(fig_dir, "initial_listing_single_pollutant_model_tier_forest.pdf")
all_png <- file.path(fig_dir, "initial_listing_all_specifications_model_tier_forest.png")
all_pdf <- file.path(fig_dir, "initial_listing_all_specifications_model_tier_forest.pdf")

ggsave(single_png, single_plot, width = 18.5, height = 12.5, dpi = 360, bg = "white")
ggsave(single_pdf, single_plot, width = 18.5, height = 12.5, device = cairo_pdf, bg = "white")
ggsave(all_png, all_spec_plot, width = 18, height = 13, dpi = 360, bg = "white")
ggsave(all_pdf, all_spec_plot, width = 18, height = 13, device = cairo_pdf, bg = "white")

write_csv(
  bind_rows(
    single_dat %>% mutate(figure = "single_pollutant_model_tier_forest"),
    all_dat %>% mutate(figure = "all_specifications_model_tier_forest")
  ) %>%
    transmute(
      figure,
      organ = organ_label,
      pollutant,
      model_tier,
      model_type,
      n,
      adverse_events,
      hazard_ratio,
      conf_low,
      conf_high,
      hr_ci
    ),
  file.path(fig_dir, "initial_listing_model_tier_forest_plot_data.csv")
)

write_csv(
  tibble(
    figure = c(
      "initial_listing_single_pollutant_model_tier_forest",
      "initial_listing_single_pollutant_model_tier_forest",
      "initial_listing_all_specifications_model_tier_forest",
      "initial_listing_all_specifications_model_tier_forest"
    ),
    path = normalizePath(c(single_png, single_pdf, all_png, all_pdf), winslash = "/", mustWork = FALSE)
  ),
  file.path(fig_dir, "initial_listing_model_tier_forest_manifest.csv")
)

message("Wrote initial-listing model tier forest plots to ", normalizePath(fig_dir, winslash = "/"))
