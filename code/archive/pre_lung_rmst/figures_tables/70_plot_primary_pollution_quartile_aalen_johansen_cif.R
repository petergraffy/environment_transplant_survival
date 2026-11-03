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
fig_dir <- file.path("output", "figures", "primary_pollution_quartile_aalen_johansen_cif")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

organ_labels <- c(HR = "Heart", KI = "Kidney", LI = "Liver", LU = "Lung")
organ_levels <- c("Heart", "Kidney", "Liver", "Lung")

pollutant_specs <- tibble(
  pollutant = c("pm25", "o3", "no2"),
  pollutant_label = c("PM[2.5]", "O[3]", "NO[2]"),
  unit_label = c("ug/m3", "ppb", "ppb"),
  exposure_col = c("pm25_waitlist_ug_m3", "o3_waitlist_ppb", "no2_waitlist"),
  event_col = c("event_pm25", "event_o3", "event_no2"),
  followup_col = c("followup_days_pm25", "followup_days_o3", "followup_days_no2"),
  end_year = c(2023L, 2023L, 2025L),
  exposure_end_date = as.Date(c("2023-12-31", "2023-12-31", "2025-12-31"))
)

quartile_labels <- c("Q1 lowest", "Q2", "Q3", "Q4 highest")
quartile_colors <- c(
  "Q1 lowest" = "#2166AC",
  "Q2" = "#67A9CF",
  "Q3" = "#F4A582",
  "Q4 highest" = "#B2182B"
)
timepoint_years <- c(1, 3, 5, 10)

log_msg <- function(...) {
  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
  flush.console()
}

pollutant_title <- function(pollutant_value, organ_value = NULL) {
  organ_text <- if (is.null(organ_value)) NULL else as.character(organ_value)
  switch(
    as.character(pollutant_value),
    "PM[2.5]" = bquote(PM[2.5]*","~.(organ_text)),
    "O[3]" = bquote(O[3]*","~.(organ_text)),
    "NO[2]" = bquote(NO[2]*","~.(organ_text)),
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
      observed_end_date,
      pollutant = spec$pollutant,
      pollutant_label = spec$pollutant_label,
      unit_label = spec$unit_label,
      exposure = .data[[spec$exposure_col]],
      followup_days = pmax(.data[[spec$followup_col]], 0.5),
      adverse = as.integer(.data[[spec$event_col]] == 1),
      competing = as.integer(
        .data[[spec$event_col]] == 0 &
          transplant_or_improvement == 1L &
          observed_end_date <= spec$exposure_end_date
      )
    ) %>%
    mutate(
      event_status = factor(
        case_when(
          adverse == 1L ~ "adverse",
          competing == 1L ~ "transplant_or_improvement",
          TRUE ~ "censor"
        ),
        levels = c("censor", "adverse", "transplant_or_improvement")
      )
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

tidy_aj <- function(dat) {
  fit <- survfit(Surv(followup_days / 365.25, event_status) ~ quartile, data = dat)
  s <- summary(fit)
  adverse_col <- match("adverse", s$states)
  transplant_col <- match("transplant_or_improvement", s$states)
  n_risk_start <- if (is.matrix(s$n.risk)) s$n.risk[, 1] else s$n.risk
  tibble(
    time = s$time,
    cif_adverse = as.numeric(s$pstate[, adverse_col]),
    cif_transplant_or_improvement = as.numeric(s$pstate[, transplant_col]),
    n_risk = as.numeric(n_risk_start),
    strata = as.character(s$strata)
  ) %>%
    mutate(
      quartile = sub("^quartile=", "", strata),
      quartile = factor(quartile, levels = quartile_labels),
      pollutant = first(dat$pollutant),
      pollutant_label = first(dat$pollutant_label),
      unit_label = first(dat$unit_label),
      organ = first(dat$organ)
    )
}

tidy_aj_times <- function(dat, times) {
  fit <- survfit(Surv(followup_days / 365.25, event_status) ~ quartile, data = dat)
  s <- summary(fit, times = times, extend = TRUE)
  adverse_col <- match("adverse", s$states)
  transplant_col <- match("transplant_or_improvement", s$states)
  n_risk_start <- if (is.matrix(s$n.risk)) s$n.risk[, 1] else s$n.risk
  tibble(
    time = s$time,
    cif_adverse = as.numeric(s$pstate[, adverse_col]),
    cif_transplant_or_improvement = as.numeric(s$pstate[, transplant_col]),
    n_risk = as.numeric(n_risk_start),
    strata = as.character(s$strata)
  ) %>%
    mutate(
      quartile = sub("^quartile=", "", strata),
      quartile = factor(quartile, levels = quartile_labels),
      pollutant = first(dat$pollutant),
      pollutant_label = first(dat$pollutant_label),
      unit_label = first(dat$unit_label),
      organ = first(dat$organ)
    )
}

log_msg("Reading primary analysis dataset")
analysis_dat <- read_csv(in_path, show_col_types = FALSE) %>%
  mutate(
    organ = factor(recode(WL_ORG, !!!organ_labels), levels = organ_levels),
    listing_year_int = as.integer(as.character(listing_year)),
    observed_end_date = as.Date(observed_end_date)
  )

plot_dat <- bind_rows(lapply(seq_len(nrow(pollutant_specs)), function(i) {
  make_pollutant_dat(analysis_dat, pollutant_specs[i, ])
}))

event_summary <- plot_dat %>%
  group_by(pollutant, pollutant_label, organ, quartile) %>%
  summarise(
    n = n(),
    adverse_events = sum(event_status == "adverse", na.rm = TRUE),
    competing_events = sum(event_status == "transplant_or_improvement", na.rm = TRUE),
    censored = sum(event_status == "censor", na.rm = TRUE),
    median_exposure = median(exposure, na.rm = TRUE),
    .groups = "drop"
  )
write_csv(event_summary, file.path(fig_dir, "pollution_quartile_aalen_johansen_event_summary.csv"))

split_dat <- split(plot_dat, list(plot_dat$pollutant, plot_dat$organ), drop = TRUE)
curve_dat <- bind_rows(lapply(split_dat, tidy_aj)) %>%
  mutate(
    pollutant_label = factor(pollutant_label, levels = pollutant_specs$pollutant_label),
    organ = factor(organ, levels = organ_levels)
  )
timepoint_dat <- bind_rows(lapply(split_dat, tidy_aj_times, times = timepoint_years)) %>%
  mutate(
    pollutant_label = factor(pollutant_label, levels = pollutant_specs$pollutant_label),
    organ = factor(organ, levels = organ_levels),
    cif_adverse_percent = 100 * cif_adverse,
    cif_transplant_or_improvement_percent = 100 * cif_transplant_or_improvement
  )

write_csv(curve_dat, file.path(fig_dir, "pollution_quartile_aalen_johansen_curve_data.csv"))
write_csv(timepoint_dat, file.path(fig_dir, "pollution_quartile_aalen_johansen_cif_at_1_3_5_10_years.csv"))

theme_cif <- function() {
  theme_minimal(base_size = 18) +
    theme(
      text = element_text(color = "grey10"),
      panel.grid = element_blank(),
      panel.border = element_rect(fill = NA, color = "grey45", linewidth = 0.55),
      axis.line = element_line(color = "black", linewidth = 0.4),
      axis.ticks = element_line(color = "black", linewidth = 0.35),
      axis.ticks.length = unit(4, "pt"),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 18),
      legend.key.width = unit(30, "pt"),
      legend.key.height = unit(13, "pt"),
      plot.title = element_text(face = "bold", size = 21),
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 16, color = "grey20"),
      panel.spacing.x = unit(12, "pt"),
      panel.spacing.y = unit(12, "pt"),
      plot.margin = margin(3, 10, 7, 8)
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
  axis_limits <- c(0, horizon)

  panel_plot <- function(pollutant_value, organ_value) {
    panel_dat <- curve_dat %>%
      filter(pollutant_label == pollutant_value, organ == organ_value, time <= horizon)

    ggplot(panel_dat, aes(x = time, y = cif_adverse, color = quartile)) +
      geom_step(linewidth = 0.86) +
      scale_color_manual(values = quartile_colors) +
      scale_x_continuous(
        breaks = breaks,
        labels = label_number(accuracy = 0.1, trim = TRUE),
        limits = axis_limits,
        expand = expansion(mult = c(0.01, 0.01))
      ) +
      scale_y_continuous(labels = label_percent(accuracy = 1), expand = expansion(mult = c(0.02, 0.04))) +
      labs(
        title = pollutant_title(pollutant_value, organ_value),
        x = "Years",
        y = "Cumulative incidence"
      ) +
      theme_cif() +
      theme(
        plot.title = element_text(face = "bold", size = 21, hjust = 0.5, margin = margin(b = 2))
      )
  }

  panel_list <- lapply(levels(curve_dat$pollutant_label), function(pollutant_value) {
    row_panels <- lapply(organ_levels, function(organ_value) {
      panel_plot(pollutant_value, organ_value)
    })
    wrap_plots(row_panels, nrow = 1)
  })

  wrap_plots(panel_list, ncol = 1) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = paste0(horizon, "-Year Cumulative Incidence of Death or Deterioration Delisting by Pollution Quartile")
    ) &
    theme(
      plot.title = element_text(face = "bold", size = 25, hjust = 0, margin = margin(b = 8)),
      plot.margin = margin(8, 14, 8, 14),
      legend.position = "bottom"
    )
}

plot_paths <- bind_rows(lapply(timepoint_years, function(horizon) {
  p <- plot_for_horizon(horizon)
  stem <- paste0("supplemental_pollution_quartile_aalen_johansen_cif_", horizon, "yr")
  png_path <- file.path(fig_dir, paste0(stem, ".png"))
  pdf_path <- file.path(fig_dir, paste0(stem, ".pdf"))
  ggsave(png_path, p, width = 22, height = 18, dpi = 360, device = "png", bg = "white")
  ggsave(pdf_path, p, width = 22, height = 18, device = grDevices::pdf, bg = "white")
  tibble(horizon_years = horizon, figure_path = c(png_path, pdf_path))
}))

file.copy(
  file.path(fig_dir, "supplemental_pollution_quartile_aalen_johansen_cif_10yr.png"),
  file.path(fig_dir, "supplemental_pollution_quartile_aalen_johansen_cif.png"),
  overwrite = TRUE
)
file.copy(
  file.path(fig_dir, "supplemental_pollution_quartile_aalen_johansen_cif_10yr.pdf"),
  file.path(fig_dir, "supplemental_pollution_quartile_aalen_johansen_cif.pdf"),
  overwrite = TRUE
)
plot_paths <- bind_rows(
  plot_paths,
  tibble(
    horizon_years = 10,
    figure_path = c(
      file.path(fig_dir, "supplemental_pollution_quartile_aalen_johansen_cif.png"),
      file.path(fig_dir, "supplemental_pollution_quartile_aalen_johansen_cif.pdf")
    )
  )
)

writeLines(
  c(
    "Figure note:",
    "Curves show nonparametric Aalen-Johansen cumulative incidence estimates for death or delisting due to deterioration, stratified by pollutant-specific waitlist-period exposure quartile.",
    "Transplant and delisting due to improvement are treated as competing events.",
    "Other waitlist exits, administrative end of follow-up, and pollutant data end are treated as censoring events."
  ),
  file.path(fig_dir, "pollution_quartile_aalen_johansen_cif_note.txt")
)

write_csv(
  plot_paths %>%
    mutate(description = paste0(horizon_years, "-year Aalen-Johansen cumulative incidence curves by waitlist-period pollutant quartile.")),
  file.path(fig_dir, "pollution_quartile_aalen_johansen_cif_manifest.csv")
)

log_msg("Wrote Aalen-Johansen CIF curves to ", normalizePath(fig_dir, winslash = "/"))
