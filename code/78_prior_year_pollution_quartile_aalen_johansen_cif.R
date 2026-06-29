#!/usr/bin/env Rscript

runtime_source <- file.path("code", "r_runtime.R")
if (!file.exists(runtime_source)) runtime_source <- "r_runtime.R"
if (file.exists(runtime_source)) {
  source(runtime_source)
  ensure_user_library()
}

suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(scales)
  library(stringr)
  library(survival)
  library(tibble)
  library(tidyr)
})

analysis_path <- file.path(
  "output",
  "primary_waitlist_period_pollution_cox",
  "primary_waitlist_period_pollution_analysis_dataset.csv.gz"
)
pollution_dir <- file.path("data", "release")
annual_pollution_dir <- file.path(pollution_dir, "air_pollution_zcta_parquet")
out_dir <- file.path("output", "prior_year_pollution_quartile_aalen_johansen_cif")
fig_dir <- file.path("output", "figures", "prior_year_pollution_quartile_aalen_johansen_cif")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

organ_labels <- c(HR = "Heart", KI = "Kidney", LI = "Liver", LU = "Lung")
organ_levels <- c("Heart", "Kidney", "Liver", "Lung")
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

clean_zip <- function(x) {
  z <- str_extract(as.character(x), "[0-9]+")
  z <- if_else(nchar(z) > 5L, substr(z, 1, 5), z)
  if_else(is.na(z), NA_character_, str_pad(z, 5, side = "left", pad = "0"))
}

parquet_files <- function(path) {
  list.files(path, pattern = "[.]parquet$", full.names = TRUE)
}

make_daily_annual <- function(path, value_col, out_col) {
  open_dataset(parquet_files(path)) %>%
    transmute(zip = zip, year = year, value = .data[[value_col]]) %>%
    group_by(zip, year) %>%
    summarise(
      n_days = sum(!is.na(value)),
      value = mean(value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(n_days >= 300, is.finite(value)) %>%
    collect() %>%
    transmute(zip = clean_zip(zip), year = as.integer(year), !!out_col := value)
}

read_prior_pollution <- function() {
  log_msg("Summarising daily PM2.5 and O3 to annual prior-year values")
  pm25 <- make_daily_annual(
    file.path(pollution_dir, "lghap_pm25_zcta_daily_parquet"),
    "pm25_ug_m3",
    "pm25_prior_ug_m3"
  )
  o3 <- make_daily_annual(
    file.path(pollution_dir, "o3_zcta_daily_parquet"),
    "o3_ppb",
    "o3_prior_ppb"
  )
  no2 <- read_parquet(file.path(annual_pollution_dir, "air_pollution_zcta_no2_annual_2005_2025.parquet")) %>%
    transmute(zip = clean_zip(zip), year = as.integer(year), no2_prior_ppb = no2)

  list(pm25 = pm25, o3 = o3, no2 = no2)
}

make_pollutant_dat <- function(dat, prior_pollution, pollutant, exposure_col, unit_label) {
  base <- dat %>%
    mutate(prior_exposure_year = listing_year_int - 1L) %>%
    left_join(prior_pollution, by = c("candidate_zip" = "zip", "prior_exposure_year" = "year")) %>%
    filter(
      is.finite(.data[[exposure_col]]),
      followup_days > 0,
      !is.na(adverse_event),
      !is.na(transplant_or_improvement)
    ) %>%
    transmute(
      organ,
      WL_ORG,
      pollutant,
      unit_label,
      exposure = .data[[exposure_col]],
      followup_years = pmax(followup_days, 0.5) / 365.25,
      event_status = factor(
        case_when(
          adverse_event == 1L ~ "adverse",
          transplant_or_improvement == 1L ~ "transplant_or_improvement",
          TRUE ~ "censor"
        ),
        levels = c("censor", "adverse", "transplant_or_improvement")
      )
    )

  cuts <- as.numeric(quantile(base$exposure, probs = seq(0, 1, 0.25), na.rm = TRUE, type = 7))
  cuts <- unique(cuts)
  if (length(cuts) < 5L) return(tibble())
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
  fit <- survfit(Surv(followup_years, event_status) ~ quartile, data = dat)
  s <- summary(fit)
  adverse_col <- match("adverse", s$states)
  transplant_col <- match("transplant_or_improvement", s$states)
  tibble(
    time = s$time,
    cif_adverse = as.numeric(s$pstate[, adverse_col]),
    cif_transplant_or_improvement = as.numeric(s$pstate[, transplant_col]),
    n_risk = as.numeric(if (is.matrix(s$n.risk)) s$n.risk[, 1] else s$n.risk),
    strata = as.character(s$strata)
  ) %>%
    mutate(
      quartile = factor(sub("^quartile=", "", strata), levels = quartile_labels),
      pollutant = first(dat$pollutant),
      unit_label = first(dat$unit_label),
      organ = first(dat$organ)
    )
}

tidy_aj_times <- function(dat, times) {
  fit <- survfit(Surv(followup_years, event_status) ~ quartile, data = dat)
  s <- summary(fit, times = times, extend = TRUE)
  adverse_col <- match("adverse", s$states)
  transplant_col <- match("transplant_or_improvement", s$states)
  tibble(
    time = s$time,
    cif_adverse = as.numeric(s$pstate[, adverse_col]),
    cif_transplant_or_improvement = as.numeric(s$pstate[, transplant_col]),
    n_risk = as.numeric(if (is.matrix(s$n.risk)) s$n.risk[, 1] else s$n.risk),
    strata = as.character(s$strata)
  ) %>%
    mutate(
      quartile = factor(sub("^quartile=", "", strata), levels = quartile_labels),
      pollutant = first(dat$pollutant),
      unit_label = first(dat$unit_label),
      organ = first(dat$organ),
      cif_adverse_percent = 100 * cif_adverse,
      cif_transplant_or_improvement_percent = 100 * cif_transplant_or_improvement
    )
}

pollutant_title <- function(pollutant_value, organ_value = NULL) {
  organ_text <- as.character(organ_value)
  switch(
    as.character(pollutant_value),
    "PM2.5" = bquote(PM[2.5]*","~.(organ_text)),
    "O3" = bquote(O[3]*","~.(organ_text)),
    "NO2" = bquote(NO[2]*","~.(organ_text)),
    as.character(pollutant_value)
  )
}

theme_cif <- function() {
  theme_minimal(base_size = 18) +
    theme(
      text = element_text(color = "grey10"),
      panel.grid = element_blank(),
      panel.border = element_rect(fill = NA, color = "grey45", linewidth = 0.55),
      axis.line = element_line(color = "black", linewidth = 0.4),
      axis.ticks = element_line(color = "black", linewidth = 0.35),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 18),
      legend.key.width = unit(30, "pt"),
      plot.title = element_text(face = "bold", size = 21),
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 16, color = "grey20"),
      panel.spacing.x = unit(12, "pt"),
      panel.spacing.y = unit(12, "pt")
    )
}

plot_for_horizon <- function(curve_dat, horizon) {
  breaks <- switch(
    as.character(horizon),
    `1` = c(0, 0.25, 0.5, 0.75, 1),
    `3` = seq(0, 3, 0.5),
    `5` = seq(0, 5, 1),
    `10` = seq(0, 10, 2),
    pretty(c(0, horizon), n = 5)
  )

  panel_plot <- function(pollutant_value, organ_value) {
    panel_dat <- curve_dat %>% filter(pollutant == pollutant_value, organ == organ_value, time <= horizon)
    ggplot(panel_dat, aes(x = time, y = cif_adverse, color = quartile)) +
      geom_step(linewidth = 0.86) +
      scale_color_manual(values = quartile_colors) +
      scale_x_continuous(
        breaks = breaks,
        labels = label_number(accuracy = 0.1, trim = TRUE),
        limits = c(0, horizon),
        expand = expansion(mult = c(0.01, 0.01))
      ) +
      scale_y_continuous(labels = label_percent(accuracy = 1), expand = expansion(mult = c(0.02, 0.04))) +
      labs(title = pollutant_title(pollutant_value, organ_value), x = "Years", y = "Cumulative incidence") +
      theme_cif() +
      theme(plot.title = element_text(hjust = 0.5, margin = margin(b = 2)))
  }

  panel_list <- lapply(c("PM2.5", "O3", "NO2"), function(pollutant_value) {
    wrap_plots(lapply(organ_levels, function(organ_value) panel_plot(pollutant_value, organ_value)), nrow = 1)
  })

  wrap_plots(panel_list, ncol = 1) +
    plot_layout(guides = "collect") +
    plot_annotation(title = paste0(horizon, "-Year Prior-Year Pollution Quartile Cumulative Incidence")) &
    theme(
      plot.title = element_text(face = "bold", size = 25, hjust = 0, margin = margin(b = 8)),
      legend.position = "bottom"
    )
}

log_msg("Reading primary deduplicated cohort")
analysis_dat <- read_csv(analysis_path, show_col_types = FALSE) %>%
  mutate(
    candidate_zip = clean_zip(candidate_zip),
    organ = factor(recode(WL_ORG, !!!organ_labels), levels = organ_levels),
    listing_year_int = as.integer(as.character(listing_year)),
    index_date = as.Date(index_date),
    observed_end_date = as.Date(observed_end_date),
    followup_days = as.numeric(observed_end_date - index_date)
  )

prior_pollution <- read_prior_pollution()
plot_dat <- bind_rows(
  make_pollutant_dat(analysis_dat, prior_pollution$pm25, "PM2.5", "pm25_prior_ug_m3", "ug/m3"),
  make_pollutant_dat(analysis_dat, prior_pollution$o3, "O3", "o3_prior_ppb", "ppb"),
  make_pollutant_dat(analysis_dat, prior_pollution$no2, "NO2", "no2_prior_ppb", "ppb")
)

event_summary <- plot_dat %>%
  group_by(pollutant, unit_label, organ, quartile) %>%
  summarise(
    n = n(),
    adverse_events = sum(event_status == "adverse"),
    competing_events = sum(event_status == "transplant_or_improvement"),
    censored = sum(event_status == "censor"),
    median_prior_exposure = median(exposure, na.rm = TRUE),
    .groups = "drop"
  )
write_csv(event_summary, file.path(out_dir, "prior_year_pollution_quartile_event_summary.csv"))

quartile_definitions <- plot_dat %>%
  group_by(pollutant, unit_label, quartile) %>%
  summarise(
    min_exposure = min(exposure, na.rm = TRUE),
    median_exposure = median(exposure, na.rm = TRUE),
    max_exposure = max(exposure, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )
write_csv(quartile_definitions, file.path(out_dir, "prior_year_pollution_national_quartile_definitions.csv"))

split_dat <- split(plot_dat, list(plot_dat$pollutant, plot_dat$organ), drop = TRUE)
curve_dat <- bind_rows(lapply(split_dat, tidy_aj)) %>%
  mutate(
    pollutant = factor(pollutant, levels = c("PM2.5", "O3", "NO2")),
    organ = factor(organ, levels = organ_levels)
  )
timepoint_dat <- bind_rows(lapply(split_dat, tidy_aj_times, times = timepoint_years)) %>%
  mutate(
    pollutant = factor(pollutant, levels = c("PM2.5", "O3", "NO2")),
    organ = factor(organ, levels = organ_levels)
  )

write_csv(curve_dat, file.path(out_dir, "prior_year_pollution_quartile_aalen_johansen_curve_data.csv"))
write_csv(timepoint_dat, file.path(out_dir, "prior_year_pollution_quartile_cif_at_1_3_5_10_years.csv"))

for (horizon in timepoint_years) {
  p <- plot_for_horizon(curve_dat, horizon)
  ggsave(
    file.path(fig_dir, paste0("prior_year_pollution_quartile_aalen_johansen_cif_", horizon, "yr.png")),
    p,
    width = 21,
    height = 18,
    dpi = 320,
    bg = "white"
  )
  ggsave(
    file.path(fig_dir, paste0("prior_year_pollution_quartile_aalen_johansen_cif_", horizon, "yr.pdf")),
    p,
    width = 21,
    height = 18,
    device = cairo_pdf,
    bg = "white"
  )
}

write_csv(
  tibble(
    figure = paste0("prior_year_pollution_quartile_aalen_johansen_cif_", timepoint_years, "yr"),
    png = normalizePath(file.path(fig_dir, paste0("prior_year_pollution_quartile_aalen_johansen_cif_", timepoint_years, "yr.png")), winslash = "/", mustWork = FALSE),
    pdf = normalizePath(file.path(fig_dir, paste0("prior_year_pollution_quartile_aalen_johansen_cif_", timepoint_years, "yr.pdf")), winslash = "/", mustWork = FALSE)
  ),
  file.path(fig_dir, "prior_year_pollution_quartile_aalen_johansen_manifest.csv")
)

log_msg("Wrote prior-year quartile Aalen-Johansen outputs to ", normalizePath(out_dir, winslash = "/"))
