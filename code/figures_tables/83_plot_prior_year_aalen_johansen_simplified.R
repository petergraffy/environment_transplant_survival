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
  library(tibble)
})

curve_path <- file.path(
  "output",
  "prior_year_pollution_quartile_aalen_johansen_cif",
  "prior_year_pollution_quartile_aalen_johansen_curve_data.csv"
)
fig_dir <- file.path("output", "figures", "prior_year_pollution_quartile_aalen_johansen_simplified")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

organ_levels <- c("Kidney", "Liver", "Heart", "Lung")
pollutant_levels_main <- c("PM2.5", "NO2")
pollutant_labels <- c(PM2.5 = "PM[2.5]", NO2 = "NO[2]", O3 = "O[3]")
quartile_labels <- c("Q1 lowest", "Q2", "Q3", "Q4 highest")
quartile_legend_labels <- c(
  "Q1 lowest" = "1st quartile (lowest)",
  "Q2" = "2nd quartile",
  "Q3" = "3rd quartile",
  "Q4 highest" = "4th quartile (highest)"
)
quartile_colors <- c(
  "Q1 lowest" = "#2166AC",
  "Q2" = "#67A9CF",
  "Q3" = "#F4A582",
  "Q4 highest" = "#B2182B"
)
horizons <- c(1, 3, 5, 10)

x_breaks_for_horizon <- function(horizon) {
  if (horizon == 1) return(c(0, 0.5, 1))
  if (horizon == 3) return(seq(0, 3, 1))
  if (horizon == 5) return(seq(0, 5, 1))
  seq(0, 10, 2)
}

x_labels_for_horizon <- function(horizon) {
  if (horizon == 1) {
    return(function(x) ifelse(x %in% c(0, 0.5, 1), c("0", "0.5", "1")[match(x, c(0, 0.5, 1))], ""))
  }
  label_number(accuracy = 1)
}

theme_aj <- function(base_size = 25) {
  theme_minimal(base_size = base_size) +
    theme(
      text = element_text(color = "grey10"),
      panel.grid = element_blank(),
      panel.border = element_rect(fill = NA, color = "black", linewidth = 0.9),
      axis.line = element_line(color = "black", linewidth = 0.55),
      axis.ticks = element_line(color = "black", linewidth = 0.5),
      axis.ticks.length = unit(5, "pt"),
      axis.text = element_text(color = "grey10", size = rel(0.9)),
      axis.title = element_text(size = rel(1.0)),
      strip.text = element_text(face = "bold", color = "grey10", size = rel(1.0), margin = margin(14, 14, 14, 14)),
      strip.text.y.right = element_text(angle = 270),
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.9),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = rel(0.9)),
      legend.key.width = unit(46, "pt"),
      panel.spacing.x = unit(16, "pt"),
      panel.spacing.y = unit(16, "pt"),
      plot.title = element_blank(),
      plot.margin = margin(12, 16, 12, 12)
    )
}

make_plot <- function(dat, horizon, title) {
  ggplot(dat %>% filter(time <= horizon), aes(x = time, y = cif_adverse, color = quartile)) +
    geom_step(linewidth = 0.95) +
    facet_grid(
      organ ~ pollutant,
      scales = "free_y",
      labeller = labeller(pollutant = as_labeller(pollutant_labels, label_parsed))
    ) +
    scale_color_manual(values = quartile_colors, breaks = quartile_labels, labels = quartile_legend_labels) +
    guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
    scale_x_continuous(
      breaks = x_breaks_for_horizon(horizon),
      limits = c(0, horizon),
      labels = x_labels_for_horizon(horizon),
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    scale_y_continuous(labels = label_percent(accuracy = 1), expand = expansion(mult = c(0.02, 0.07))) +
    labs(x = "Years after listing", y = "Cumulative incidence") +
    theme_aj() +
    theme(
      strip.switch.pad.grid = unit(0, "pt")
    )
}

curve_dat <- read_csv(curve_path, show_col_types = FALSE) %>%
  mutate(
    pollutant = factor(pollutant, levels = c("PM2.5", "O3", "NO2")),
    organ = factor(organ, levels = organ_levels),
    quartile = factor(quartile, levels = quartile_labels)
  ) %>%
  filter(time <= max(horizons))

main_dat <- curve_dat %>%
  filter(pollutant %in% pollutant_levels_main) %>%
  mutate(pollutant = factor(as.character(pollutant), levels = pollutant_levels_main))

o3_dat <- curve_dat %>%
  filter(pollutant == "O3") %>%
  mutate(pollutant = factor(as.character(pollutant), levels = "O3"))

manifest <- list()

for (horizon in horizons) {
  horizon_label <- paste0(horizon, "yr")

  main_plot <- make_plot(
    main_dat,
    horizon,
    NULL
  )
  o3_plot <- make_plot(
    o3_dat,
    horizon,
    NULL
  )

  main_png <- file.path(fig_dir, paste0("prior_year_pm25_no2_aalen_johansen_abcd_", horizon_label, ".png"))
  main_pdf <- file.path(fig_dir, paste0("prior_year_pm25_no2_aalen_johansen_abcd_", horizon_label, ".pdf"))
  o3_png <- file.path(fig_dir, paste0("supplemental_prior_year_o3_aalen_johansen_", horizon_label, ".png"))
  o3_pdf <- file.path(fig_dir, paste0("supplemental_prior_year_o3_aalen_johansen_", horizon_label, ".pdf"))

  ggsave(main_png, main_plot, width = 16, height = 22, dpi = 360, bg = "white")
  ggsave(main_pdf, main_plot, width = 16, height = 22, device = cairo_pdf, bg = "white")
  ggsave(o3_png, o3_plot, width = 10.8, height = 22, dpi = 360, bg = "white")
  ggsave(o3_pdf, o3_plot, width = 10.8, height = 22, device = cairo_pdf, bg = "white")

  manifest[[length(manifest) + 1]] <- tibble(
    figure = c(
      paste0("prior_year_pm25_no2_aalen_johansen_abcd_", horizon_label),
      paste0("supplemental_prior_year_o3_aalen_johansen_", horizon_label)
    ),
    horizon_years = horizon,
    png = normalizePath(c(main_png, o3_png), winslash = "/", mustWork = FALSE),
    pdf = normalizePath(c(main_pdf, o3_pdf), winslash = "/", mustWork = FALSE)
  )
}

write_csv(
  bind_rows(
    main_dat %>% mutate(display_figure = "main_pm25_no2_abcd"),
    o3_dat %>% mutate(display_figure = "supplement_o3")
  ),
  file.path(fig_dir, "prior_year_aalen_johansen_simplified_curve_data.csv")
)
write_csv(
  bind_rows(manifest),
  file.path(fig_dir, "prior_year_aalen_johansen_simplified_manifest.csv")
)

message("Wrote simplified prior-year Aalen-Johansen figures to ", normalizePath(fig_dir, winslash = "/"))
