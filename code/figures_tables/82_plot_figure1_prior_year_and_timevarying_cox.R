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
  library(patchwork)
  library(readr)
  library(scales)
  library(stringr)
  library(tibble)
})

baseline_path <- file.path("output", "prior_year_pollution_cox_svi", "prior_year_pollution_cox_svi_results.csv")
timevarying_path <- file.path("output", "timevarying_pollution_cox_svi", "timevarying_pollution_cox_svi_results.csv")
fig_dir <- file.path("output", "figures", "figure1_cox_baseline_timevarying")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

pollutant_levels <- c("pm25", "o3", "no2")
pollutant_labels_plotmath <- c(
  pm25 = 'PM[2.5]~"per 5 ug/m"^3',
  o3 = 'O[3]~"per 10 ppb"',
  no2 = 'NO[2]~"per 10 ppb"'
)
pollutant_colors <- c(pm25 = "#0072B2", o3 = "#009E73", no2 = "#D55E00")
organ_levels <- c("Kidney", "Liver", "Heart", "Lung")
model_levels <- c("1 year before listing", "Waitlist-period exposure")

format_ci <- function(hr, low, high) {
  sprintf("%.2f (%.2f-%.2f)", hr, low, high)
}

theme_figure1 <- function(base_size = 17) {
  theme_minimal(base_size = base_size) +
    theme(
      text = element_text(color = "grey10"),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(linewidth = 0.28, color = "grey86"),
      axis.title.y = element_blank(),
      axis.text = element_text(color = "grey15", size = rel(0.95)),
      axis.title.x = element_text(size = rel(1.0), margin = margin(t = 10)),
      strip.text = element_text(face = "bold", color = "grey10", size = rel(1.0)),
      strip.background = element_rect(fill = "grey94", color = "grey68", linewidth = 0.6),
      panel.border = element_rect(fill = NA, color = "grey72", linewidth = 0.55),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = rel(0.9)),
      legend.key.width = unit(24, "pt"),
      plot.title = element_blank(),
      plot.margin = margin(8, 12, 8, 8)
    )
}

baseline <- read_csv(baseline_path, show_col_types = FALSE) %>%
  transmute(
    model = "1 year before listing",
    organ = organ_label,
    pollutant = exposure,
    hazard_ratio,
    conf_low,
    conf_high,
    n,
    adverse_events
  )

timevarying <- read_csv(timevarying_path, show_col_types = FALSE) %>%
  transmute(
    model = "Waitlist-period exposure",
    organ = organ_label,
    pollutant,
    hazard_ratio,
    conf_low,
    conf_high,
    n = candidate_organ_episodes,
    adverse_events
  )

plot_dat <- bind_rows(baseline, timevarying) %>%
  mutate(
    model = factor(model, levels = model_levels),
    organ = factor(organ, levels = rev(organ_levels)),
    pollutant = factor(pollutant, levels = pollutant_levels),
    hr_ci = format_ci(hazard_ratio, conf_low, conf_high),
    label_left = hazard_ratio > 1.75,
    label_x = if_else(label_left, conf_low / 1.04, conf_high * 1.04),
    label_hjust = if_else(label_left, 1, 0)
  )

p <- ggplot(plot_dat, aes(x = hazard_ratio, y = organ, color = pollutant, group = pollutant)) +
  geom_vline(xintercept = 1, linewidth = 0.4, color = "grey35") +
  geom_errorbar(
    aes(xmin = conf_low, xmax = conf_high),
    orientation = "y",
    width = 0.18,
    linewidth = 0.62,
    position = position_dodge(width = 0.72)
  ) +
  geom_point(size = 3.0, position = position_dodge(width = 0.72)) +
  geom_text(
    aes(x = label_x, label = hr_ci, hjust = label_hjust),
    position = position_dodge(width = 0.72),
    size = 4.0,
    color = "grey10",
    show.legend = FALSE
  ) +
  facet_wrap(~model, nrow = 1) +
  scale_color_manual(
    values = pollutant_colors,
    breaks = pollutant_levels,
    labels = parse(text = pollutant_labels_plotmath[pollutant_levels])
  ) +
  scale_x_log10(
    breaks = c(0.8, 1.0, 1.25, 1.5, 2.0),
    labels = label_number(accuracy = 0.01),
    limits = c(0.78, 2.18),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  coord_cartesian(clip = "off") +
  labs(
    x = "Cause-specific hazard ratio, log scale"
  ) +
  theme_figure1()

png_path <- file.path(fig_dir, "figure1_prior_year_and_timevarying_cox.png")
pdf_path <- file.path(fig_dir, "figure1_prior_year_and_timevarying_cox.pdf")
ggsave(png_path, p, width = 15.5, height = 7.2, dpi = 360, bg = "white")
ggsave(pdf_path, p, width = 15.5, height = 7.2, device = cairo_pdf, bg = "white")

write_csv(
  plot_dat %>%
    transmute(model, organ, pollutant, n, adverse_events, hazard_ratio, conf_low, conf_high, hr_ci),
  file.path(fig_dir, "figure1_prior_year_and_timevarying_cox_plot_table.csv")
)
write_csv(
  tibble(
    figure = "figure1_prior_year_and_timevarying_cox",
    path = normalizePath(c(png_path, pdf_path), winslash = "/", mustWork = FALSE)
  ),
  file.path(fig_dir, "figure1_prior_year_and_timevarying_cox_manifest.csv")
)

message("Wrote Figure 1 Cox comparison to ", normalizePath(fig_dir, winslash = "/"))
