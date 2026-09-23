#!/usr/bin/env Rscript

runtime_source <- file.path("code", "r_runtime.R")
if (!file.exists(runtime_source)) runtime_source <- "r_runtime.R"
if (file.exists(runtime_source)) {
  source(runtime_source)
  ensure_user_library()
}

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(tidyr)
})

out_dir <- file.path("output", "tables", "updated_methods_sensitivity_models")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

organ_levels <- c("Heart", "Kidney", "Liver", "Lung")
pollutant_levels <- c("PM2.5", "O3", "NO2")
model_columns <- c(
  "Primary model",
  "No listing center strata",
  "Single-organ candidates only",
  "Primary model + multi-organ status",
  "Multipollutant PM2.5 + NO2",
  "Multipollutant PM2.5 + NO2 + O3"
)

format_hr <- function(hr, low, high) {
  if_else(
    is.na(hr) | is.na(low) | is.na(high),
    "",
    sprintf("%.2f (%.2f-%.2f)", hr, low, high)
  )
}

pollutant_from_text <- function(x) {
  case_when(
    grepl("pm25|PM2[.]5", x, ignore.case = TRUE) ~ "PM2.5",
    grepl("o3|O3", x, ignore.case = TRUE) ~ "O3",
    grepl("no2|NO2", x, ignore.case = TRUE) ~ "NO2",
    TRUE ~ NA_character_
  )
}

make_wide_table <- function(primary, sensitivity, analysis_name, primary_label) {
  if (!"exposure" %in% names(primary)) primary$exposure <- NA_character_
  if (!"pollutant" %in% names(primary)) primary$pollutant <- NA_character_
  if (!"exposure_label" %in% names(primary)) primary$exposure_label <- NA_character_

  primary_long <- primary %>%
    transmute(
      analysis = analysis_name,
      model = primary_label,
      organ_label,
      pollutant = pollutant_from_text(coalesce(exposure, pollutant, exposure_label)),
      result = format_hr(hazard_ratio, conf_low, conf_high)
    )

  sensitivity_long <- sensitivity %>%
    filter(analysis == analysis_name) %>%
    mutate(
      pollutant = pollutant_from_text(coalesce(exposure_label, exposure)),
      model = recode(
        model,
        "Single-pollutant model without listing center" = "No listing center strata",
        "Single-pollutant model restricted to single-organ candidates" = "Single-organ candidates only",
        "Single-pollutant model additionally adjusted for multi-organ candidate status" = "Primary model + multi-organ status",
        "Multipollutant model: PM2.5 + NO2" = "Multipollutant PM2.5 + NO2",
        "Multipollutant model: PM2.5 + NO2 + O3" = "Multipollutant PM2.5 + NO2 + O3"
      )
    ) %>%
    transmute(
      analysis,
      model,
      organ_label,
      pollutant,
      result = format_hr(hazard_ratio, conf_low, conf_high)
    )

  bind_rows(primary_long, sensitivity_long) %>%
    mutate(
      organ_label = factor(organ_label, levels = organ_levels),
      pollutant = factor(pollutant, levels = pollutant_levels),
      model = factor(
        model,
        levels = model_columns
      )
    ) %>%
    filter(!is.na(organ_label), !is.na(pollutant), !is.na(model)) %>%
    distinct(organ_label, pollutant, model, .keep_all = TRUE) %>%
    select(organ_label, pollutant, model, result) %>%
    pivot_wider(names_from = model, values_from = result, values_fill = "") %>%
    arrange(organ_label, pollutant) %>%
    mutate(
      Organ = as.character(organ_label),
      Pollutant = as.character(pollutant),
      Organ = if_else(duplicated(Organ), "", Organ)
    ) %>%
    select(Organ, Pollutant, all_of(model_columns), -organ_label, -pollutant)
}

baseline_primary <- read_csv(
  file.path("output", "prior_year_pollution_cox_svi", "prior_year_pollution_cox_svi_results.csv"),
  show_col_types = FALSE
) %>%
  mutate(pollutant = exposure)

timevarying_primary <- read_csv(
  file.path("output", "timevarying_pollution_cox_svi", "timevarying_pollution_cox_svi_results.csv"),
  show_col_types = FALSE
) %>%
  rename(n = candidate_organ_episodes)

sensitivity <- read_csv(
  file.path("output", "updated_methods_sensitivity_models", "updated_methods_sensitivity_results.csv"),
  show_col_types = FALSE
)

baseline_table <- make_wide_table(
  baseline_primary,
  sensitivity,
  "Baseline cause-specific Cox",
  "Primary model"
)

timevarying_table <- make_wide_table(
  timevarying_primary,
  sensitivity,
  "Time-varying cause-specific Cox",
  "Primary model"
)

combined_table <- bind_rows(
  baseline_table %>% mutate(Analysis = "Baseline cause-specific Cox", .before = 1),
  timevarying_table %>% mutate(Analysis = "Time-varying cause-specific Cox", .before = 1)
)

write_csv(baseline_table, file.path(out_dir, "supplemental_baseline_sensitivity_wide_table.csv"))
write_tsv(baseline_table, file.path(out_dir, "supplemental_baseline_sensitivity_wide_table.tsv"))
write_csv(timevarying_table, file.path(out_dir, "supplemental_timevarying_sensitivity_wide_table.csv"))
write_tsv(timevarying_table, file.path(out_dir, "supplemental_timevarying_sensitivity_wide_table.tsv"))
write_csv(combined_table, file.path(out_dir, "supplemental_sensitivity_wide_table_set.csv"))
write_tsv(combined_table, file.path(out_dir, "supplemental_sensitivity_wide_table_set.tsv"))

notes <- c(
  "Supplemental wide sensitivity table set.",
  "Primary model rows are the updated single-pollutant primary models for each analysis: 1-year exposure before listing for baseline Cox models and monthly waitlist-period exposure for time-varying Cox models.",
  "The PM2.5 + NO2 multipollutant model does not include O3; therefore O3 cells are blank in that column.",
  "Values are hazard ratios with 95% confidence intervals."
)
writeLines(notes, file.path(out_dir, "supplemental_sensitivity_wide_table_set_notes.txt"))

message("Wrote updated wide sensitivity table set to ", normalizePath(out_dir, winslash = "/"))
