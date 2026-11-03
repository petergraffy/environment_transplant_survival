#!/usr/bin/env Rscript

runtime_source <- file.path("code", "r_runtime.R")
if (!file.exists(runtime_source)) runtime_source <- "r_runtime.R"
if (file.exists(runtime_source)) {
  source(runtime_source)
  ensure_user_library()
}

suppressPackageStartupMessages({
  library(broom)
  library(dplyr)
  library(readr)
  library(survival)
  library(tibble)
})

in_path <- file.path(
  "output",
  "primary_waitlist_period_pollution_cox",
  "primary_waitlist_period_pollution_analysis_dataset.csv.gz"
)
analysis_start_year <- as.integer(Sys.getenv("FG_ANALYSIS_START_YEAR", "2015"))
analysis_label <- paste0(analysis_start_year, "plus")
out_dir <- file.path("output", paste0("primary_waitlist_period_pollution_finegray_", analysis_label))
model_result_dir <- file.path(out_dir, "model_results")
log_dir <- file.path("logs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(model_result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

target_organs <- c("HR", "KI", "LI", "LU")
organ_labels <- c(HR = "Heart", KI = "Kidney", LI = "Liver", LU = "Lung")

exposure_specs <- list(
  pm25 = list(
    term = "pm25_waitlist_5ug",
    label = "PM2.5 waitlist-period mean per 5 ug/m3",
    followup = "followup_days_pm25",
    end_date = as.Date("2023-12-31"),
    end_year = 2023L
  ),
  o3 = list(
    term = "o3_waitlist_10ppb",
    label = "Ozone waitlist-period mean per 10 ppb",
    followup = "followup_days_o3",
    end_date = as.Date("2023-12-31"),
    end_year = 2023L
  ),
  no2 = list(
    term = "no2_waitlist_10unit",
    label = "NO2 waitlist-period mean per 10-unit increment",
    followup = "followup_days_no2",
    end_date = as.Date("2025-12-31"),
    end_year = 2025L
  )
)

progress_path <- file.path(out_dir, "finegray_progress.csv")

log_msg <- function(...) {
  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
  flush.console()
}

append_progress <- function(row) {
  row <- as_tibble(row)
  write_csv(row, progress_path, append = file.exists(progress_path))
}

safe_nrow <- function(x) {
  if (is.null(x)) return(NA_real_)
  nr <- tryCatch(nrow(x), error = function(e) NA_integer_)
  if (is.null(nr) || length(nr) == 0L || is.na(nr)) return(NA_real_)
  as.numeric(nr)
}

safe_error <- function(x) {
  if (is.null(x) || length(x) == 0L) return(NA_character_)
  as.character(x[[1]])
}

make_progress <- function(total, label) {
  pb <- txtProgressBar(min = 0, max = total, style = 3)
  attr(pb, "label") <- label
  pb
}

tick_progress <- function(pb, value, detail = "") {
  setTxtProgressBar(pb, value)
  if (nzchar(detail)) message("\n", attr(pb, "label"), " ", value, "/", attr(pb, "total"), " | ", detail)
  invisible(NULL)
}

format_status <- function(dat, spec) {
  event_date <- dat$observed_end_date
  within_exposure_followup <- !is.na(event_date) & event_date <= spec$end_date
  case_when(
    within_exposure_followup & dat$event_type == "death_or_deterioration_delist" ~ "adverse",
    within_exposure_followup & dat$event_type == "transplant_or_improvement" ~ "transplant_or_improvement",
    TRUE ~ "censor"
  )
}

fit_finegray <- function(dat, org, exposure_name, spec) {
  result_path <- file.path(model_result_dir, paste0(tolower(org), "_", exposure_name, "_finegray.csv"))
  if (file.exists(result_path) && Sys.getenv("FG_OVERWRITE", "FALSE") != "TRUE") {
    log_msg("Skipping existing Fine-Gray result ", org, " | ", exposure_name)
    return(result_path)
  }

  vars_needed <- c(
    spec$followup, "event_type", "observed_end_date", spec$term,
    "age", "sex", "race", "listing_year", "organ_score", "listing_center"
  )
  model_dat <- dat %>%
    filter(WL_ORG == org, complete.cases(across(all_of(vars_needed))), .data[[spec$followup]] > 0) %>%
    mutate(
      .fg_time = pmax(.data[[spec$followup]], 0.5),
      .exposure = .data[[spec$term]],
      .fg_status = factor(
        format_status(cur_data_all(), spec),
        levels = c("censor", "adverse", "transplant_or_improvement")
      ),
      sex = factor(sex),
      race = factor(race),
      listing_year = factor(listing_year),
      listing_center = factor(listing_center)
    ) %>%
    droplevels()

  source_n <- nrow(model_dat)
  source_adverse <- sum(model_dat$.fg_status == "adverse")
  source_competing <- sum(model_dat$.fg_status == "transplant_or_improvement")
  if (source_n == 0L || source_adverse < 50L || source_competing < 50L) {
    write_csv(tibble(), result_path)
    return(result_path)
  }

  log_msg(
    "Fine-Gray full cohort start ", org, " | ", exposure_name,
    " n=", source_n,
    " adverse=", source_adverse,
    " competing=", source_competing
  )
  started_at <- Sys.time()
  append_progress(tibble(
    organ = org,
    exposure = exposure_name,
    status = "started",
    started_at = as.character(started_at),
    finished_at = NA_character_,
    elapsed_seconds = NA_real_,
    n = source_n,
    adverse_events = source_adverse,
    competing_events = source_competing,
    expanded_rows = NA_real_,
    error = NA_character_
  ))

  fg_dat <- NULL
  fit <- NULL
  fit_error <- NULL
  elapsed <- system.time({
    fg_dat <- tryCatch(
      finegray(
        Surv(.fg_time, .fg_status) ~
          age + sex + race + listing_year + organ_score + listing_center + .exposure,
        data = model_dat,
        etype = "adverse"
      ),
      error = function(e) {
        fit_error <<- conditionMessage(e)
        NULL
      }
    )
    if (!is.null(fg_dat)) {
      form <- as.formula(paste(
        "Surv(fgstart, fgstop, fgstatus) ~",
        paste(c(".exposure", "age", "sex", "race", "listing_year", "organ_score", "strata(listing_center)"), collapse = " + ")
      ))
      fit <- tryCatch(
        coxph(form, data = fg_dat, weights = fgwt, ties = "efron", x = FALSE, y = FALSE),
        error = function(e) {
          fit_error <<- conditionMessage(e)
          NULL
        }
      )
    }
  })

  finished_at <- Sys.time()
  if (is.null(fit)) {
    result <- tibble(
      organ = org,
      organ_label = recode(org, !!!organ_labels),
      exposure = exposure_name,
      exposure_label = spec$label,
      exposure_window = "observed_waitlist_period_mean",
      exposure_data_end_year = spec$end_year,
      endpoint = "death_or_deterioration_delist_subdistribution",
      n = source_n,
      adverse_events = source_adverse,
      competing_transplant_or_improvement = source_competing,
      expanded_rows = safe_nrow(fg_dat),
      subdistribution_hazard_ratio = NA_real_,
      conf_low = NA_real_,
      conf_high = NA_real_,
      p_value = NA_real_,
      elapsed_seconds = unname(elapsed[["elapsed"]]),
      adjustment_set = "age + sex + race + listing_year + baseline_organ_score + strata(listing_center)",
      variance_estimator = "survival_finegray_weighted_cox_model_based",
      center_adjustment = "stratified_baseline_hazard_by_listing_center",
      sampled = FALSE,
      error = safe_error(fit_error)
    )
    write_csv(result, result_path)
  } else {
    result <- tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
      filter(term == ".exposure") %>%
      transmute(
        organ = org,
        organ_label = recode(org, !!!organ_labels),
        exposure = exposure_name,
        exposure_label = spec$label,
        exposure_window = "observed_waitlist_period_mean",
        exposure_data_end_year = spec$end_year,
        endpoint = "death_or_deterioration_delist_subdistribution",
        n = source_n,
        adverse_events = source_adverse,
        competing_transplant_or_improvement = source_competing,
        expanded_rows = safe_nrow(fg_dat),
        subdistribution_hazard_ratio = estimate,
        conf_low = conf.low,
        conf_high = conf.high,
        p_value = p.value,
        elapsed_seconds = unname(elapsed[["elapsed"]]),
        adjustment_set = "age + sex + race + listing_year + baseline_organ_score + strata(listing_center)",
        variance_estimator = "survival_finegray_weighted_cox_model_based",
        center_adjustment = "stratified_baseline_hazard_by_listing_center",
        sampled = FALSE,
        error = NA_character_
      )
    write_csv(result, result_path)
  }

  append_progress(tibble(
    organ = org,
    exposure = exposure_name,
    status = if (is.null(fit)) "error" else "completed",
    started_at = as.character(started_at),
    finished_at = as.character(finished_at),
    elapsed_seconds = unname(elapsed[["elapsed"]]),
      n = source_n,
      adverse_events = source_adverse,
      competing_events = source_competing,
      expanded_rows = safe_nrow(fg_dat),
      error = safe_error(fit_error)
    ))

  log_msg(
    "Fine-Gray full cohort done ", org, " | ", exposure_name,
    " elapsed_seconds=", round(unname(elapsed[["elapsed"]]), 1),
    " expanded_rows=", safe_nrow(fg_dat)
  )
  result_path
}

log_msg("Reading primary waitlist-period pollution analysis dataset")
dat <- read_csv(in_path, show_col_types = FALSE) %>%
  mutate(
    observed_end_date = as.Date(observed_end_date),
    listing_year_integer = as.integer(as.character(listing_year)),
    pm25_waitlist_5ug = pm25_waitlist_ug_m3 / 5,
    o3_waitlist_10ppb = o3_waitlist_ppb / 10,
    no2_waitlist_10unit = no2_waitlist / 10
  ) %>%
  filter(listing_year_integer >= analysis_start_year)

log_msg("Fine-Gray cohort restricted to listing year >= ", analysis_start_year, "; rows=", nrow(dat))

jobs <- expand.grid(organ = target_organs, exposure = names(exposure_specs), stringsAsFactors = FALSE)
pb <- make_progress(nrow(jobs), "Fine-Gray models")
attr(pb, "total") <- nrow(jobs)
result_paths <- character()
for (i in seq_len(nrow(jobs))) {
  org <- jobs$organ[[i]]
  exposure_name <- jobs$exposure[[i]]
  tick_progress(pb, i - 1L, paste(org, exposure_name, "starting"))
  result_paths <- c(result_paths, fit_finegray(dat, org, exposure_name, exposure_specs[[exposure_name]]))
  tick_progress(pb, i, paste(org, exposure_name, "finished"))
  invisible(gc())
}
close(pb)

finegray_results <- bind_rows(lapply(result_paths[file.exists(result_paths)], read_csv, show_col_types = FALSE))
write_csv(finegray_results, file.path(out_dir, "primary_waitlist_period_pollution_finegray_results.csv"))
log_msg("Wrote full-cohort primary Fine-Gray outputs to ", out_dir)
