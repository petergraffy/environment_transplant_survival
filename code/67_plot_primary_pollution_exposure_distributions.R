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
  library(tibble)
  library(tidyr)
})

in_path <- file.path(
  "output",
  "primary_waitlist_period_pollution_cox",
  "primary_waitlist_period_pollution_analysis_dataset.csv.gz"
)
fig_dir <- file.path("output", "figures", "primary_pollution_exposure_distributions")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

organ_labels <- c(HR = "Heart", KI = "Kidney", LI = "Liver", LU = "Lung")
organ_levels <- c("Heart", "Kidney", "Liver", "Lung")

pollutant_specs <- tibble(
  pollutant = c("pm25", "o3", "no2"),
  pollutant_label = c("PM[2.5]", "O[3]", "NO[2]"),
  unit_label = c("\u00b5g/m\u00b3", "ppb", "ppb"),
  exposure_col = c("pm25_waitlist_ug_m3", "o3_waitlist_ppb", "no2_waitlist"),
  event_col = c("event_pm25", "event_o3", "event_no2"),
  followup_col = c("followup_days_pm25", "followup_days_o3", "followup_days_no2"),
  end_year = c(2023L, 2023L, 2025L)
)

log_msg <- function(...) {
  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
  flush.console()
}

log_msg("Reading primary analysis dataset")
analysis_dat <- read_csv(in_path, show_col_types = FALSE) %>%
  mutate(
    organ = factor(recode(WL_ORG, !!!organ_labels), levels = organ_levels),
    listing_year_int = as.integer(as.character(listing_year))
  )

plot_dat <- bind_rows(lapply(seq_len(nrow(pollutant_specs)), function(i) {
  spec <- pollutant_specs[i, ]
  analysis_dat %>%
    filter(
      listing_year_int <= spec$end_year,
      .data[[spec$followup_col]] > 0,
      is.finite(.data[[spec$exposure_col]]),
      !is.na(.data[[spec$event_col]])
    ) %>%
    transmute(
      organ,
      pollutant = spec$pollutant,
      pollutant_label = spec$pollutant_label,
      unit_label = spec$unit_label,
      exposure = .data[[spec$exposure_col]],
      outcome = if_else(.data[[spec$event_col]] == 1, "Death/deterioration delisting", "Other/censored")
    )
}))

plot_dat <- plot_dat %>%
  mutate(
    pollutant_label = factor(pollutant_label, levels = pollutant_specs$pollutant_label),
    outcome = factor(outcome, levels = c("Other/censored", "Death/deterioration delisting"))
  )

summary_dat <- plot_dat %>%
  group_by(pollutant, pollutant_label, unit_label, organ, outcome) %>%
  summarise(
    n = n(),
    mean = mean(exposure, na.rm = TRUE),
    sd = sd(exposure, na.rm = TRUE),
    median = median(exposure, na.rm = TRUE),
    p25 = as.numeric(quantile(exposure, 0.25, na.rm = TRUE)),
    p75 = as.numeric(quantile(exposure, 0.75, na.rm = TRUE)),
    .groups = "drop"
  )
write_csv(summary_dat, file.path(fig_dir, "primary_pollution_exposure_distribution_summary.csv"))

plot_limits <- plot_dat %>%
  group_by(pollutant_label) %>%
  summarise(
    lo = as.numeric(quantile(exposure, 0.01, na.rm = TRUE)),
    hi = as.numeric(quantile(exposure, 0.99, na.rm = TRUE)),
    .groups = "drop"
  )

plot_dat <- plot_dat %>%
  left_join(plot_limits, by = "pollutant_label") %>%
  mutate(exposure_plot = pmin(pmax(exposure, lo), hi))

theme_exposure <- function() {
  theme_minimal(base_size = 12) +
    theme(
      text = element_text(color = "grey10"),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(linewidth = 0.25, color = "grey88"),
      axis.title.x = element_text(margin = margin(t = 8)),
      axis.title.y = element_text(margin = margin(r = 8)),
      strip.text = element_text(face = "bold"),
      strip.background = element_rect(fill = "grey94", color = "grey72", linewidth = 0.4),
      legend.position = "top",
      legend.title = element_blank(),
      plot.title = element_text(face = "bold", size = 16),
      plot.margin = margin(8, 10, 8, 8)
    )
}

median_dat <- plot_dat %>%
  group_by(pollutant_label, organ, outcome) %>%
  summarise(median_exposure = median(exposure, na.rm = TRUE), .groups = "drop")

p <- ggplot(plot_dat, aes(x = exposure_plot, fill = outcome, color = outcome)) +
  geom_density(
    alpha = 0.32,
    linewidth = 0.75,
    adjust = 1.05
  ) +
  geom_vline(
    data = median_dat,
    aes(xintercept = median_exposure, color = outcome),
    linewidth = 0.45,
    linetype = "dashed",
    show.legend = FALSE
  ) +
  facet_grid(pollutant_label ~ organ, scales = "free_y", labeller = labeller(pollutant_label = label_parsed)) +
  scale_fill_manual(values = c("Other/censored" = "#9ECAE1", "Death/deterioration delisting" = "#F1695B")) +
  scale_color_manual(values = c("Other/censored" = "#3182BD", "Death/deterioration delisting" = "#CB181D")) +
  scale_x_continuous(labels = label_number(accuracy = 0.1)) +
  labs(
    title = "Waitlist-Period Pollution Exposure Distributions",
    x = "Waitlist-period average exposure",
    y = "Density"
  ) +
  theme_exposure()

png_path <- file.path(fig_dir, "supplemental_waitlist_pollution_exposure_distributions.png")
pdf_path <- file.path(fig_dir, "supplemental_waitlist_pollution_exposure_distributions.pdf")

ggsave(png_path, p, width = 11.5, height = 9.2, dpi = 500, bg = "white")
ggsave(pdf_path, p, width = 11.5, height = 9.2, device = grDevices::pdf, bg = "white")

write_csv(
  tibble(
    figure_path = c(png_path, pdf_path),
    description = "Supplemental boxplots of waitlist-period pollutant exposure by organ and pollutant-specific adverse outcome status."
  ),
  file.path(fig_dir, "primary_pollution_exposure_distribution_manifest.csv")
)

log_msg("Wrote supplemental exposure distribution figure to ", normalizePath(fig_dir, winslash = "/"))
