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
  library(survival)
  library(tibble)
  library(tidyr)
})

tv_results_path <- file.path("output", "formal_waitlist_environment_timevarying", "timevarying_cox_exposure_results.csv")
formal_dataset_path <- file.path("output", "formal_waitlist_environment", "formal_waitlist_environment_analysis_dataset.csv.gz")
out_dir <- file.path("output", "figures", "formal_waitlist_environment_timevarying")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

max_km_years <- 5
max_km_days <- max_km_years * 365.25

organ_levels <- c("Heart", "Kidney", "Liver", "Lung")
exposure_levels <- c(
  "Interval Tmax per 5 C",
  "Interval maximum relative humidity per 10 pct",
  "Interval PM2.5 per 5 ug/m3",
  "Interval ozone per 10 ppb",
  "Interval NO2 per 10 units"
)
exposure_labels <- c(
  "Interval Tmax per 5 C" = "Tmax",
  "Interval maximum relative humidity per 10 pct" = "Max RH",
  "Interval PM2.5 per 5 ug/m3" = "PM2.5",
  "Interval ozone per 10 ppb" = "Ozone",
  "Interval NO2 per 10 units" = "NO2"
)

theme_pub <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom"
    )
}

log_msg <- function(...) {
  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
  flush.console()
}

log_msg("Plotting time-varying Cox forest")
tv <- read_csv(tv_results_path, show_col_types = FALSE) %>%
  mutate(
    organ_label = factor(organ_label, levels = organ_levels),
    exposure_short = recode(exposure, !!!exposure_labels),
    exposure_short = factor(exposure_short, levels = rev(c("Tmax", "Max RH", "PM2.5", "Ozone", "NO2")))
  )

forest_path <- file.path(out_dir, "timevarying_cox_exposure_forest.png")
forest <- ggplot(tv, aes(x = hazard_ratio, y = exposure_short)) +
  geom_vline(xintercept = 1, color = "grey45", linewidth = 0.35) +
  geom_errorbarh(aes(xmin = conf_low, xmax = conf_high), height = 0.18, linewidth = 0.45) +
  geom_point(size = 1.9) +
  scale_x_log10() +
  facet_wrap(~organ_label, ncol = 2) +
  labs(
    title = "Time-Varying Cox Models of Environment and Adverse Waitlist Outcome",
    subtitle = "Start-stop status-history intervals; center-stratified baseline hazards and robust patient clustering",
    x = "Hazard ratio for death/deterioration delisting",
    y = NULL
  ) +
  theme_pub()
ggsave(forest_path, forest, width = 11, height = 7, dpi = 240, bg = "white")

log_msg("Building Kaplan-Meier curves by exposure quartile")
km_specs <- tribble(
  ~exposure, ~term, ~label,
  "tmax", "tmax_c_waitlist_mean_exact_daily", "Tmax",
  "pm25", "pm25_waitlist_ug_m3", "PM2.5",
  "o3", "o3_waitlist_ppb", "Ozone",
  "no2", "no2_waitlist", "NO2"
)

make_quartile <- function(x) {
  qs <- quantile(x, probs = seq(0, 1, 0.25), na.rm = TRUE, names = FALSE)
  qs <- unique(qs)
  if (length(qs) < 3L) return(factor(rep(NA_character_, length(x))))
  cut(x, breaks = qs, include.lowest = TRUE, labels = paste0("Q", seq_len(length(qs) - 1L)))
}

tidy_survfit <- function(fit, org, exposure_name, exposure_label) {
  s <- summary(fit)
  strata <- if (is.null(s$strata)) rep("All", length(s$time)) else as.character(s$strata)
  tibble(
    organ = org,
    organ_label = recode(org, HR = "Heart", KI = "Kidney", LI = "Liver", LU = "Lung"),
    exposure = exposure_name,
    exposure_label = exposure_label,
    quartile = sub("^exposure_quartile=", "", strata),
    time_days = s$time,
    survival = s$surv,
    cumulative_adverse = 1 - s$surv,
    n_risk = s$n.risk,
    n_event = s$n.event
  )
}

base <- read_csv(formal_dataset_path, show_col_types = FALSE) %>%
  filter(WL_ORG %in% c("HR", "KI", "LI", "LU"), followup_days > 0) %>%
  mutate(adverse_event = as.integer(adverse_event == 1L))

km_curves <- bind_rows(lapply(seq_len(nrow(km_specs)), function(i) {
  spec <- km_specs[i, ]
  bind_rows(lapply(c("HR", "KI", "LI", "LU"), function(org) {
    dat <- base %>%
      filter(WL_ORG == org, !is.na(.data[[spec$term]])) %>%
      mutate(exposure_quartile = make_quartile(.data[[spec$term]])) %>%
      filter(!is.na(exposure_quartile))
    if (nrow(dat) == 0L || n_distinct(dat$exposure_quartile) < 2L) return(tibble())
    log_msg("KM ", org, " | ", spec$label, " n=", nrow(dat), " adverse=", sum(dat$adverse_event))
    fit <- survfit(Surv(followup_days, adverse_event) ~ exposure_quartile, data = dat)
    tidy_survfit(fit, org, spec$exposure, spec$label)
  }))
}))

write_csv(km_curves, file.path("output", "formal_waitlist_environment_timevarying", "km_by_exposure_quartile_curves.csv"))

km_plot_dat <- km_curves %>%
  filter(time_days <= max_km_days) %>%
  mutate(
    years = time_days / 365.25,
    organ_label = factor(organ_label, levels = organ_levels),
    quartile = factor(quartile, levels = paste0("Q", 1:4))
  )

km_survival_path <- file.path(out_dir, "kaplan_meier_survival_by_exposure_quartile.png")
km_surv <- ggplot(km_plot_dat, aes(x = years, y = survival, color = quartile)) +
  geom_step(linewidth = 0.7) +
  facet_grid(exposure_label ~ organ_label) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  scale_color_brewer(palette = "Dark2", na.translate = FALSE) +
  labs(
    title = "Kaplan-Meier Adverse-Event-Free Survival by Exposure Quartile",
    subtitle = paste0(
      "Restricted to ", max_km_years,
      " years; death/deterioration delisting is the event; transplant, improvement, other removals, and administrative end are censored"
    ),
    x = "Years from waitlist start",
    y = "Adverse-event-free survival",
    color = "Quartile"
  ) +
  theme_pub()
ggsave(km_survival_path, km_surv, width = 14, height = 9, dpi = 240, bg = "white")

km_adverse_path <- file.path(out_dir, "kaplan_meier_cumulative_adverse_by_exposure_quartile.png")
km_adv <- ggplot(km_plot_dat, aes(x = years, y = cumulative_adverse, color = quartile)) +
  geom_step(linewidth = 0.7) +
  facet_grid(exposure_label ~ organ_label) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_color_brewer(palette = "Dark2", na.translate = FALSE) +
  labs(
    title = "Kaplan-Meier Cumulative Adverse Event by Exposure Quartile",
    subtitle = paste0("Restricted to ", max_km_years, " years; 1 - KM survival; competing exits are treated as censoring"),
    x = "Years from waitlist start",
    y = "1 - KM survival",
    color = "Quartile"
  ) +
  theme_pub()
ggsave(km_adverse_path, km_adv, width = 14, height = 9, dpi = 240, bg = "white")

write_csv(
  tibble(
    figure = c("timevarying_cox_exposure_forest", "kaplan_meier_survival_by_exposure_quartile", "kaplan_meier_cumulative_adverse_by_exposure_quartile"),
    path = c(forest_path, km_survival_path, km_adverse_path)
  ),
  file.path(out_dir, "timevarying_cox_km_figure_manifest.csv")
)

log_msg("Wrote time-varying Cox and KM figures to ", out_dir)
