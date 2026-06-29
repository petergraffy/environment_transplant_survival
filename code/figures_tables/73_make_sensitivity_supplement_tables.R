#!/usr/bin/env Rscript

runtime_source <- file.path("code", "r_runtime.R")
if (!file.exists(runtime_source)) runtime_source <- "r_runtime.R"
if (file.exists(runtime_source)) {
  source(runtime_source)
  ensure_user_library()
}

suppressPackageStartupMessages({
  library(dplyr)
  library(gt)
  library(readr)
  library(stringr)
  library(tibble)
  library(tidyr)
})

primary_path <- file.path(
  "output",
  "primary_waitlist_period_pollution_cox",
  "primary_waitlist_period_pollution_cox_results.csv"
)
sensitivity_path <- file.path(
  "output",
  "primary_waitlist_period_pollution_cox_sensitivity",
  "primary_waitlist_period_pollution_cox_sensitivity_results.csv"
)
out_dir <- file.path("output", "tables", "primary_waitlist_period_pollution_sensitivity")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

pollutant_labels <- c(
  pm25 = "PM2.5",
  o3 = "O3",
  no2 = "NO2"
)

model_labels <- c(
  primary = "Primary model",
  no_organ_score = "No organ score",
  no_listing_center = "No listing center strata",
  primary_plus_svi = "Primary model + SVI proxy",
  single_organ_candidates = "Single-organ candidates only",
  primary_plus_multi_organ = "Primary model + multi-organ status",
  prior_1y_no_score_no_center = "Prior-year exposure, no organ score or center strata",
  multipollutant_pm25_no2 = "Multipollutant PM2.5 + NO2",
  multipollutant_pm25_no2_o3 = "Multipollutant PM2.5 + NO2 + O3"
)

model_order <- names(model_labels)
organ_order <- c("Heart", "Kidney", "Liver", "Lung")
pollutant_order <- unname(pollutant_labels)

fmt_int <- function(x) {
  formatC(as.integer(round(x)), format = "d", big.mark = ",")
}

fmt_hr_ci <- function(hr, low, high) {
  ifelse(
    is.na(hr) | is.na(low) | is.na(high),
    "",
    paste0(sprintf("%.2f", hr), " (", sprintf("%.2f", low), "-", sprintf("%.2f", high), ")")
  )
}

fmt_n_events <- function(n, events) {
  ifelse(
    is.na(n) | is.na(events),
    "",
    paste0(fmt_int(n), " / ", fmt_int(events))
  )
}

primary_results <- read_csv(primary_path, show_col_types = FALSE) %>%
  mutate(
    sensitivity = "primary",
    sensitivity_label = model_labels[["primary"]]
  )

sensitivity_results <- read_csv(sensitivity_path, show_col_types = FALSE)

all_results <- bind_rows(primary_results, sensitivity_results) %>%
  mutate(
    model = factor(sensitivity, levels = model_order),
    model_label = factor(recode(sensitivity, !!!model_labels), levels = unname(model_labels)),
    organ_label = factor(organ_label, levels = organ_order),
    pollutant = factor(recode(exposure, !!!pollutant_labels), levels = pollutant_order),
    hr_ci = fmt_hr_ci(hazard_ratio, conf_low, conf_high),
    n_events = fmt_n_events(n, adverse_events)
  ) %>%
  filter(!is.na(model), !is.na(organ_label), !is.na(pollutant))

hr_table <- all_results %>%
  select(Organ = organ_label, Pollutant = pollutant, model_label, hr_ci) %>%
  pivot_wider(names_from = model_label, values_from = hr_ci) %>%
  mutate(across(where(is.character), ~ replace_na(.x, ""))) %>%
  arrange(Organ, Pollutant)

n_event_table <- all_results %>%
  select(Organ = organ_label, Pollutant = pollutant, model_label, n_events) %>%
  pivot_wider(names_from = model_label, values_from = n_events) %>%
  mutate(across(where(is.character), ~ replace_na(.x, ""))) %>%
  arrange(Organ, Pollutant)

model_spec_table <- all_results %>%
  mutate(adjustment_group = if_else(as.character(organ_label) == "Kidney", "Kidney", "Heart, liver, and lung")) %>%
  distinct(
    sensitivity,
    model_label,
    adjustment_group,
    exposure_window,
    adjustment_set,
    center_adjustment,
    variance_estimator
  ) %>%
  mutate(
    exposure_window = as.character(exposure_window),
    center_adjustment = as.character(center_adjustment),
    exposure_window = case_when(
      sensitivity == "prior_1y_no_score_no_center" ~ "Mean exposure during the calendar year before listing",
      sensitivity == "single_organ_candidates" ~ "Day-weighted mean exposure during observed waitlist follow-up; restricted to candidates listed for one organ group",
      str_detect(sensitivity, "^multipollutant_") ~ "Day-weighted mean exposure during observed waitlist follow-up; common data horizon through 2023",
      exposure_window == "observed_waitlist_period_mean" ~ "Day-weighted mean exposure during observed waitlist follow-up",
      exposure_window == "prior_year_mean" ~ "Mean exposure during the calendar year before listing",
      exposure_window == "one year prior listing" ~ "Mean exposure during the calendar year before listing",
      TRUE ~ str_replace_all(exposure_window, "_", " ")
    ),
    center_adjustment = case_when(
      center_adjustment == "stratified_baseline_hazard_by_listing_center" ~ "Listing center strata",
      center_adjustment == "not_included" ~ "Not included",
      center_adjustment == "none" ~ "Not included",
      TRUE ~ str_replace_all(center_adjustment, "_", " ")
    ),
    variance_estimator = case_when(
      variance_estimator == "model_based" ~ "Model-based",
      TRUE ~ str_replace_all(variance_estimator, "_", " ")
    )
  ) %>%
  arrange(factor(sensitivity, levels = model_order)) %>%
  transmute(
    Model = as.character(model_label),
    `Applicable organs` = adjustment_group,
    `Exposure window` = exposure_window,
    `Adjustment set` = adjustment_set,
    `Center adjustment` = center_adjustment,
    `Variance estimator` = variance_estimator
  )

write_csv(hr_table, file.path(out_dir, "supplemental_sensitivity_hr_table.csv"))
write_csv(n_event_table, file.path(out_dir, "supplemental_sensitivity_sample_size_events_table.csv"))
write_csv(model_spec_table, file.path(out_dir, "supplemental_sensitivity_model_specifications.csv"))

hr_gt <- hr_table %>%
  gt(rowname_col = "Pollutant", groupname_col = "Organ") %>%
  tab_header(title = "Sensitivity Analyses for Cause-Specific Cox Models") %>%
  tab_source_note(source_note = "Values are hazard ratios (95% CIs) for death or delisting due to clinical deterioration.") %>%
  tab_source_note(source_note = "PM2.5 is scaled per 5-ug/m3 increase; O3 and NO2 are scaled per 10-ppb increase.") %>%
  tab_options(
    table.font.names = "Arial",
    table.font.size = px(12),
    heading.title.font.size = px(16),
    row_group.font.weight = "bold",
    column_labels.font.weight = "bold",
    data_row.padding = px(4)
  )

n_event_gt <- n_event_table %>%
  gt(rowname_col = "Pollutant", groupname_col = "Organ") %>%
  tab_header(title = "Analytic Sample Size and Adverse Events in Sensitivity Analyses") %>%
  tab_source_note(source_note = "Values are person waitlist episodes / adverse events.") %>%
  tab_options(
    table.font.names = "Arial",
    table.font.size = px(12),
    heading.title.font.size = px(16),
    row_group.font.weight = "bold",
    column_labels.font.weight = "bold",
    data_row.padding = px(4)
  )

model_spec_gt <- model_spec_table %>%
  gt() %>%
  tab_header(title = "Sensitivity Model Specifications") %>%
  tab_options(
    table.font.names = "Arial",
    table.font.size = px(12),
    heading.title.font.size = px(16),
    column_labels.font.weight = "bold",
    data_row.padding = px(4)
  )

gtsave(hr_gt, file.path(out_dir, "supplemental_sensitivity_hr_table.html"))
gtsave(n_event_gt, file.path(out_dir, "supplemental_sensitivity_sample_size_events_table.html"))
gtsave(model_spec_gt, file.path(out_dir, "supplemental_sensitivity_model_specifications.html"))

writeLines(
  c(
    "Supplemental sensitivity table notes:",
    "The primary model used day-weighted waitlist-period pollutant exposure and adjusted for age, sex, race, listing year, organ-specific baseline severity, and listing-center strata. For kidney candidates, baseline severity adjustment was represented by no dialysis time, dialysis duration, and diabetes rather than a composite organ score.",
    "Sensitivity models removed organ-specific severity adjustment, removed listing-center strata, added an ACS-derived ZCTA-level SVI proxy, restricted the cohort to single-organ candidates, added multi-organ candidate status, or used prior-year exposure without organ score or center strata.",
    "PM2.5 is scaled per 5-ug/m3 increase; O3 and NO2 are scaled per 10-ppb increase."
  ),
  file.path(out_dir, "supplemental_sensitivity_table_notes.txt")
)

message("Wrote sensitivity supplement tables to ", normalizePath(out_dir, winslash = "/"))
