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
  library(survival)
  library(tibble)
  library(tidyr)
})

in_path <- file.path(
  "output",
  "primary_waitlist_period_pollution_cox",
  "primary_waitlist_period_pollution_analysis_dataset.csv.gz"
)
fig_dir <- file.path("output", "figures", "primary_pollution_quartile_km_curves")
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

quartile_labels <- c("Q1 lowest", "Q2", "Q3", "Q4 highest")
quartile_colors <- c(
  "Q1 lowest" = "#2166AC",
  "Q2" = "#67A9CF",
  "Q3" = "#F4A582",
  "Q4 highest" = "#B2182B"
)
timepoint_years <- c(1, 3, 5, 10)

format_cut <- function(x) {
  format(round(x, 2), nsmall = 2, trim = TRUE)
}

pollutant_title <- function(pollutant_value, organ_value) {
  switch(
    as.character(pollutant_value),
    "PM[2.5]" = bquote(PM[2.5]~"/ "~.(as.character(organ_value))),
    "O[3]" = bquote(O[3]~"/ "~.(as.character(organ_value))),
    "NO[2]" = bquote(NO[2]~"/ "~.(as.character(organ_value))),
    as.character(pollutant_value)
  )
}

make_pollutant_dat <- function(dat, spec) {
  base <- dat %>%
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
      followup_days = pmax(.data[[spec$followup_col]], 0.5),
      event = as.integer(.data[[spec$event_col]] == 1)
    )

  cuts <- as.numeric(quantile(base$exposure, probs = seq(0, 1, 0.25), na.rm = TRUE, type = 7))
  cuts[1] <- -Inf
  cuts[length(cuts)] <- Inf

  base %>%
    mutate(
      quartile = cut(
        exposure,
        breaks = cuts,
        labels = quartile_labels,
        include.lowest = TRUE,
        right = TRUE
      ),
      quartile = factor(quartile, levels = quartile_labels)
    )
}

tidy_survfit <- function(fit) {
  s <- summary(fit)
  tibble(
    time = s$time,
    survival = s$surv,
    std_error = s$std.err,
    n_risk = s$n.risk,
    n_event = s$n.event,
    n_censor = s$n.censor,
    strata = s$strata
  )
}

tidy_survfit_times <- function(fit, times) {
  s <- summary(fit, times = times, extend = TRUE)
  tibble(
    time = s$time,
    survival = s$surv,
    std_error = s$std.err,
    n_risk = s$n.risk,
    n_event = s$n.event,
    n_censor = s$n.censor,
    strata = s$strata
  )
}

log_msg("Reading primary analysis dataset")
analysis_dat <- read_csv(in_path, show_col_types = FALSE) %>%
  mutate(
    organ = factor(recode(WL_ORG, !!!organ_labels), levels = organ_levels),
    listing_year_int = as.integer(as.character(listing_year))
  )

plot_dat <- bind_rows(lapply(seq_len(nrow(pollutant_specs)), function(i) {
  make_pollutant_dat(analysis_dat, pollutant_specs[i, ])
}))

quartile_summary <- plot_dat %>%
  group_by(pollutant, pollutant_label, unit_label, quartile) %>%
  summarise(
    n = n(),
    min = min(exposure, na.rm = TRUE),
    median = median(exposure, na.rm = TRUE),
    max = max(exposure, na.rm = TRUE),
    .groups = "drop"
  )
write_csv(quartile_summary, file.path(fig_dir, "pollution_quartile_definitions.csv"))

curve_dat <- bind_rows(lapply(split(plot_dat, list(plot_dat$pollutant, plot_dat$organ), drop = TRUE), function(dat) {
  fit <- survfit(Surv(followup_days / 365.25, event) ~ quartile, data = dat)
  tidy_survfit(fit) %>%
    mutate(
      quartile = sub("^quartile=", "", strata),
      quartile = factor(quartile, levels = quartile_labels),
      pollutant = first(dat$pollutant),
      pollutant_label = first(dat$pollutant_label),
      unit_label = first(dat$unit_label),
      organ = first(dat$organ),
      adverse_probability = 1 - survival
    )
}))

timepoint_dat <- bind_rows(lapply(split(plot_dat, list(plot_dat$pollutant, plot_dat$organ), drop = TRUE), function(dat) {
  fit <- survfit(Surv(followup_days / 365.25, event) ~ quartile, data = dat)
  tidy_survfit_times(fit, timepoint_years) %>%
    mutate(
      quartile = sub("^quartile=", "", strata),
      quartile = factor(quartile, levels = quartile_labels),
      pollutant = first(dat$pollutant),
      pollutant_label = first(dat$pollutant_label),
      unit_label = first(dat$unit_label),
      organ = first(dat$organ),
      adverse_probability = 1 - survival
    )
}))

curve_dat <- curve_dat %>%
  mutate(
    pollutant_label = factor(pollutant_label, levels = pollutant_specs$pollutant_label),
    organ = factor(organ, levels = organ_levels)
  )

timepoint_dat <- timepoint_dat %>%
  mutate(
    pollutant_label = factor(pollutant_label, levels = pollutant_specs$pollutant_label),
    organ = factor(organ, levels = organ_levels),
    survival_percent = survival * 100,
    adverse_probability_percent = adverse_probability * 100
  )

risk_table_times <- c(0, 0.25, 0.5, 0.75, 1, 1.5, 2, 2.5, 3, 4, 5, 6, 7, 8, 9, 10)
risk_table_dat <- tidyr::crossing(
  plot_dat %>% distinct(pollutant, pollutant_label, organ, quartile),
  time = risk_table_times
) %>%
  left_join(plot_dat, by = c("pollutant", "pollutant_label", "organ", "quartile"), relationship = "many-to-many") %>%
  group_by(pollutant, pollutant_label, organ, quartile, time) %>%
  summarise(
    at_risk = sum(followup_days / 365.25 >= time, na.rm = TRUE),
    censored = sum(followup_days / 365.25 <= time & event == 0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    pollutant_label = factor(pollutant_label, levels = pollutant_specs$pollutant_label),
    organ = factor(organ, levels = organ_levels),
    quartile = factor(quartile, levels = quartile_labels),
    label = paste0(at_risk, "\n", censored)
  )
write_csv(risk_table_dat, file.path(fig_dir, "pollution_quartile_km_risk_censor_table.csv"))

event_summary <- plot_dat %>%
  group_by(pollutant, pollutant_label, unit_label, organ, quartile) %>%
  summarise(
    n = n(),
    adverse_events = sum(event, na.rm = TRUE),
    median_exposure = median(exposure, na.rm = TRUE),
    .groups = "drop"
  )
write_csv(event_summary, file.path(fig_dir, "pollution_quartile_km_event_summary.csv"))
write_csv(curve_dat, file.path(fig_dir, "pollution_quartile_km_curve_data.csv"))
write_csv(timepoint_dat, file.path(fig_dir, "pollution_quartile_km_survival_at_1_3_5_10_years.csv"))

theme_km <- function() {
  theme_minimal(base_size = 16) +
    theme(
      text = element_text(color = "grey10"),
      panel.grid = element_blank(),
      panel.border = element_rect(fill = NA, color = "grey45", linewidth = 0.55),
      axis.line = element_line(color = "black", linewidth = 0.4),
      axis.ticks = element_line(color = "black", linewidth = 0.35),
      axis.ticks.length = unit(4, "pt"),
      strip.text = element_text(face = "bold", size = 15),
      strip.background = element_rect(fill = "grey94", color = "grey45", linewidth = 0.55),
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(size = 14),
      plot.title = element_text(face = "bold", size = 21),
      axis.title = element_text(size = 16),
      panel.spacing.x = unit(14, "pt"),
      panel.spacing.y = unit(10, "pt"),
      plot.margin = margin(8, 10, 8, 8)
    )
}

theme_risk_table <- function() {
  theme_minimal(base_size = 15) +
    theme(
      text = element_text(color = "grey10"),
      panel.grid = element_blank(),
      panel.border = element_rect(fill = NA, color = "grey45", linewidth = 0.55),
      axis.line.x = element_line(color = "black", linewidth = 0.35),
      axis.ticks.x = element_line(color = "black", linewidth = 0.3),
      axis.ticks.y = element_blank(),
      axis.title = element_blank(),
      axis.text.x = element_text(size = 13.2, color = "grey15"),
      axis.text.y = element_text(size = 13.2, color = "grey15", margin = margin(r = 6)),
      strip.text = element_blank(),
      strip.background = element_blank(),
      plot.margin = margin(0, 8, 4, 8)
    )
}

plot_for_horizon <- function(horizon) {
  breaks <- switch(
    as.character(horizon),
    `1` = c(0, 0.25, 0.5, 0.75, 1),
    `3` = seq(0, 3, 0.5),
    `5` = seq(0, 5, 1),
    `10` = seq(0, 10, 2),
    pretty(c(0, horizon), n = 5)
  )
  axis_limits <- c(-0.075 * horizon, 1.075 * horizon)

  panel_plot <- function(pollutant_value, organ_value) {
    curve_panel_dat <- curve_dat %>%
      filter(pollutant_label == pollutant_value, organ == organ_value, time <= horizon)
    table_panel_dat <- risk_table_dat %>%
      filter(pollutant_label == pollutant_value, organ == organ_value, time %in% breaks, time <= horizon) %>%
      mutate(quartile = factor(quartile, levels = rev(quartile_labels)))

    curve_panel <- ggplot(curve_panel_dat, aes(x = time, y = survival, color = quartile)) +
      geom_step(linewidth = 0.82) +
      scale_color_manual(values = quartile_colors) +
      scale_x_continuous(
        breaks = breaks,
        labels = NULL,
        limits = axis_limits,
        expand = expansion(mult = c(0, 0))
      ) +
      scale_y_continuous(labels = label_percent(accuracy = 1), expand = expansion(mult = c(0.04, 0.02))) +
      labs(
        title = pollutant_title(pollutant_value, organ_value),
        y = "KM survival"
      ) +
      theme_km() +
      theme(
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        plot.title = element_text(face = "bold", size = 17, hjust = 0.5, margin = margin(b = 4)),
        plot.margin = margin(4, 4, 1, 4)
      )

    table_panel <- ggplot(table_panel_dat, aes(x = time, y = quartile, label = label, color = quartile)) +
      geom_text(size = 3.9, lineheight = 0.9, show.legend = FALSE) +
      scale_color_manual(values = quartile_colors) +
      scale_x_continuous(
        breaks = breaks,
        labels = label_number(accuracy = 0.1, trim = TRUE),
        limits = axis_limits,
        expand = expansion(mult = c(0, 0))
      ) +
      scale_y_discrete(drop = FALSE) +
      coord_cartesian(clip = "off") +
      labs(x = "Years", y = NULL) +
      theme_risk_table()

    curve_panel / table_panel + plot_layout(heights = c(2.9, 1.35))
  }

  panel_list <- lapply(levels(curve_dat$pollutant_label), function(pollutant_value) {
    row_panels <- lapply(organ_levels, function(organ_value) {
      panel_plot(pollutant_value, organ_value)
    })
    wrap_plots(row_panels, nrow = 1)
  })

  wrap_plots(panel_list, ncol = 1) +
    plot_annotation(
      title = paste0(horizon, "-Year Adverse-Event-Free Waitlist Survival by Pollution Quartile"),
      caption = "Risk table cells show at risk over censored"
    ) &
    theme(
      plot.title = element_text(face = "bold", size = 24, hjust = 0, margin = margin(b = 8)),
      plot.caption = element_text(hjust = 0, size = 14, color = "grey20", margin = margin(t = 6)),
      plot.margin = margin(8, 14, 8, 14)
    )
}

plot_paths <- bind_rows(lapply(timepoint_years, function(horizon) {
  p <- plot_for_horizon(horizon)
  stem <- paste0("supplemental_pollution_quartile_km_curves_", horizon, "yr")
  png_path <- file.path(fig_dir, paste0(stem, ".png"))
  pdf_path <- file.path(fig_dir, paste0(stem, ".pdf"))
  ggsave(png_path, p, width = 24, height = 30, dpi = 320, device = "png", bg = "white")
  ggsave(pdf_path, p, width = 24, height = 30, device = grDevices::pdf, bg = "white")
  tibble(horizon_years = horizon, figure_path = c(png_path, pdf_path))
}))

# Preserve the 10-year file name used by earlier drafts.
file.copy(
  file.path(fig_dir, "supplemental_pollution_quartile_km_curves_10yr.png"),
  file.path(fig_dir, "supplemental_pollution_quartile_km_curves.png"),
  overwrite = TRUE
)
file.copy(
  file.path(fig_dir, "supplemental_pollution_quartile_km_curves_10yr.pdf"),
  file.path(fig_dir, "supplemental_pollution_quartile_km_curves.pdf"),
  overwrite = TRUE
)
plot_paths <- bind_rows(
  plot_paths,
  tibble(
    horizon_years = 10,
    figure_path = c(
      file.path(fig_dir, "supplemental_pollution_quartile_km_curves.png"),
      file.path(fig_dir, "supplemental_pollution_quartile_km_curves.pdf")
    )
  )
)

writeLines(
  c(
    "Figure note:",
    "Curves show unadjusted Kaplan-Meier survival estimates through 10 years for remaining free of death or delisting due to deterioration, stratified by pollutant-specific waitlist-period exposure quartile.",
    "Transplant, delisting due to improvement, other exits, administrative end of follow-up, and pollutant data end are treated as censoring events to match the cause-specific Cox model.",
    "These are cause-specific Kaplan-Meier curves, not competing-risk cumulative incidence functions."
  ),
  file.path(fig_dir, "pollution_quartile_km_curves_note.txt")
)

write_csv(
  plot_paths %>%
    mutate(description = paste0(horizon_years, "-year cause-specific Kaplan-Meier survival curves by waitlist-period pollutant quartile.")),
  file.path(fig_dir, "pollution_quartile_km_curve_manifest.csv")
)

log_msg("Wrote pollution quartile KM curves to ", normalizePath(fig_dir, winslash = "/"))
