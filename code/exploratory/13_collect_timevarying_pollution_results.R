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
  library(stringr)
})

out_dir <- file.path("output", "organ_specific_adverse_waitlist")
step_dirs <- file.path(out_dir, c(
  "timevarying_steps_full_seq",
  "timevarying_steps_li_lu_seq",
  "timevarying_steps_ki_seq",
  "timevarying_steps_curve_hr",
  "timevarying_steps_curve_li_lu",
  "timevarying_steps_curve_ki"
))
step_dirs <- step_dirs[dir.exists(step_dirs)]

read_tagged <- function(paths, suffix) {
  paths <- paths[str_detect(basename(paths), suffix)]
  if (length(paths) == 0L) return(tibble())
  bind_rows(lapply(paths, function(path) {
    read_csv(path, col_types = cols(.default = col_character()), show_col_types = FALSE) %>%
      mutate(step_file = basename(path), step_dir = basename(dirname(path)), .before = 1)
  })) %>%
    type_convert()
}

all_files <- unlist(lapply(step_dirs, function(path) list.files(path, full.names = TRUE, pattern = "[.]csv$")), use.names = FALSE)

linear <- read_tagged(all_files, "_linear[.]csv$")
contrast <- read_tagged(all_files, "_contrast[.]csv$")
curve <- read_tagged(all_files, "_curve[.]csv$")
progress <- bind_rows(lapply(file.path(step_dirs, "progress.csv"), function(path) {
  if (!file.exists(path)) return(tibble())
  read_csv(path, col_types = cols(.default = col_character()), show_col_types = FALSE) %>%
    mutate(step_dir = basename(dirname(path)), .before = 1)
})) %>%
  type_convert()

write_csv(linear, file.path(out_dir, "organ_specific_timevarying_pollution_linear_cox_results_all_seq.csv"))
write_csv(contrast, file.path(out_dir, "organ_specific_timevarying_pollution_spline_iqr_contrasts_all_seq.csv"))
write_csv(curve, file.path(out_dir, "organ_specific_timevarying_pollution_spline_curves_all_seq.csv"))
write_csv(progress, file.path(out_dir, "organ_specific_timevarying_pollution_progress_all_seq.csv"))

if (nrow(linear) > 0L) {
  linear_summary <- linear %>%
    filter(exposure_form == "linear", term %in% c("pm25_waitlist_cum_5ug", "o3_waitlist_cum_10ppb", "no2_waitlist_cum_10unit", "pm25_waitlist_5ug", "o3_waitlist_10ppb", "no2_waitlist_10unit")) %>%
    transmute(
      organ,
      exposure_window,
      model_family = if_else(str_detect(model, "^combined"), "combined", "single"),
      pollutant,
      term,
      hazard_ratio = estimate,
      conf_low = conf.low,
      conf_high = conf.high,
      p_value = p.value,
      n_intervals,
      adverse_events,
      concordance
    ) %>%
    arrange(organ, exposure_window, pollutant, model_family)
  write_csv(linear_summary, file.path(out_dir, "organ_specific_timevarying_pollution_linear_summary_all_seq.csv"))
}

if (nrow(contrast) > 0L) {
  contrast_summary <- contrast %>%
    filter(contrast == "p75_vs_p25") %>%
    mutate(model_family = if_else(str_detect(model, "^combined"), "combined", "single")) %>%
    select(
      organ, exposure_window, exposure_form, model_family, pollutant, p25, p75,
      hazard_ratio, conf_low, conf_high, p_value, wald_chisq, df,
      n_intervals, adverse_events, concordance
    ) %>%
    arrange(organ, exposure_window, pollutant, model_family, exposure_form)
  write_csv(contrast_summary, file.path(out_dir, "organ_specific_timevarying_pollution_spline_iqr_summary_all_seq.csv"))
}

message("Combined ", length(all_files), " step files from ", length(step_dirs), " directories.")
