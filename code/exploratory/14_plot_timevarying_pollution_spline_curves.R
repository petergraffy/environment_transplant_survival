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
})

out_dir <- file.path("output", "organ_specific_adverse_waitlist")
figure_dir <- file.path("output", "figures", "timevarying_pollution_splines")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

curve_path <- file.path(out_dir, "organ_specific_timevarying_pollution_spline_curves_all_seq.csv")
if (!file.exists(curve_path)) {
  stop("Missing spline curve file: ", curve_path, call. = FALSE)
}

curves <- read_csv(curve_path, show_col_types = FALSE) %>%
  mutate(
    organ_label = recode(organ, HR = "Heart", KI = "Kidney", LI = "Liver", LU = "Lung"),
    pollutant_label = recode(
      pollutant,
      "PM2.5" = "PM2.5 (per 5 ug/m3)",
      "Ozone" = "Ozone (per 10 ppb)",
      "NO2" = "NO2 (per 10-unit increase)"
    ),
    model_family = if_else(str_detect(model, "^combined"), "Combined pollutants", "Single pollutant"),
    window_label = recode(
      exposure_window,
      cumulative = "Cumulative waitlist exposure",
      lagged_annual = "Lagged annual waitlist exposure"
    ),
    organ_label = factor(organ_label, levels = c("Heart", "Kidney", "Liver", "Lung")),
    pollutant_label = factor(
      pollutant_label,
      levels = c("PM2.5 (per 5 ug/m3)", "Ozone (per 10 ppb)", "NO2 (per 10-unit increase)")
    ),
    model_family = factor(model_family, levels = c("Single pollutant", "Combined pollutants"))
  ) %>%
  filter(is.finite(exposure_value), is.finite(hazard_ratio))

if (nrow(curves) == 0L) stop("Spline curve file has no usable rows.", call. = FALSE)

clip_curve <- function(dat) {
  dat %>%
    mutate(
      conf_low_plot = pmax(conf_low, 0.25, na.rm = TRUE),
      conf_high_plot = pmin(conf_high, 4, na.rm = TRUE),
      hazard_ratio_plot = pmin(pmax(hazard_ratio, 0.25), 4)
    )
}

plot_window <- function(window_value) {
  dat <- curves %>%
    filter(exposure_window == window_value) %>%
    clip_curve()
  window_title <- unique(dat$window_label)
  p <- ggplot(dat, aes(x = exposure_value, y = hazard_ratio_plot, color = model_family, fill = model_family)) +
    geom_hline(yintercept = 1, linewidth = 0.35, color = "gray45") +
    geom_ribbon(aes(ymin = conf_low_plot, ymax = conf_high_plot), alpha = 0.12, color = NA) +
    geom_line(linewidth = 0.75) +
    facet_grid(organ_label ~ pollutant_label, scales = "free_x") +
    scale_y_log10(
      breaks = c(0.5, 0.75, 1, 1.5, 2, 3, 4),
      labels = label_number(accuracy = 0.01),
      limits = c(0.25, 4)
    ) +
    scale_color_manual(values = c("Single pollutant" = "#1f6f8b", "Combined pollutants" = "#b6403a")) +
    scale_fill_manual(values = c("Single pollutant" = "#1f6f8b", "Combined pollutants" = "#b6403a")) +
    labs(
      title = paste0(window_title, ": adjusted spline curves for death/delistment"),
      subtitle = "Hazard ratios are relative to the median observed exposure in each organ/model; ribbons are pointwise 95% intervals clipped to the plotted range.",
      x = "Exposure on model scale",
      y = "Adjusted hazard ratio",
      color = NULL,
      fill = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      legend.position = "top",
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(size = 8)
    )

  path <- file.path(figure_dir, paste0("timevarying_spline_curves_", window_value, "_single_vs_combined.png"))
  ggsave(path, p, width = 13.5, height = 9, dpi = 320)
  path
}

plot_pollutant <- function(pollutant_value) {
  dat <- curves %>%
    filter(pollutant == pollutant_value) %>%
    clip_curve()
  safe_pollutant <- str_to_lower(str_replace_all(pollutant_value, "[^A-Za-z0-9]+", ""))
  p <- ggplot(dat, aes(x = exposure_value, y = hazard_ratio_plot, color = model_family, fill = model_family)) +
    geom_hline(yintercept = 1, linewidth = 0.35, color = "gray45") +
    geom_ribbon(aes(ymin = conf_low_plot, ymax = conf_high_plot), alpha = 0.12, color = NA) +
    geom_line(linewidth = 0.75) +
    facet_grid(organ_label ~ window_label, scales = "free_x") +
    scale_y_log10(
      breaks = c(0.5, 0.75, 1, 1.5, 2, 3, 4),
      labels = label_number(accuracy = 0.01),
      limits = c(0.25, 4)
    ) +
    scale_color_manual(values = c("Single pollutant" = "#1f6f8b", "Combined pollutants" = "#b6403a")) +
    scale_fill_manual(values = c("Single pollutant" = "#1f6f8b", "Combined pollutants" = "#b6403a")) +
    labs(
      title = paste0(unique(dat$pollutant_label), ": adjusted time-varying spline curves"),
      subtitle = "Hazard ratios are relative to the median observed exposure in each organ/model; ribbons are pointwise 95% intervals clipped to the plotted range.",
      x = "Exposure on model scale",
      y = "Adjusted hazard ratio",
      color = NULL,
      fill = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      legend.position = "top",
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold")
    )

  path <- file.path(figure_dir, paste0("timevarying_spline_curves_", safe_pollutant, ".png"))
  ggsave(path, p, width = 10.5, height = 8, dpi = 320)
  path
}

paths <- c(
  plot_window("cumulative"),
  plot_window("lagged_annual"),
  unlist(lapply(c("PM2.5", "Ozone", "NO2"), plot_pollutant), use.names = FALSE)
)

write_csv(tibble(figure_path = paths), file.path(figure_dir, "timevarying_spline_curve_figures.csv"))
message("Wrote ", length(paths), " spline curve figures to ", figure_dir)
