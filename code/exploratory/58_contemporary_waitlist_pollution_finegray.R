#!/usr/bin/env Rscript

runtime_source <- file.path("code", "r_runtime.R")
if (!file.exists(runtime_source)) runtime_source <- "r_runtime.R"
if (file.exists(runtime_source)) {
  source(runtime_source)
  ensure_user_library()
}

suppressPackageStartupMessages({
    library(arrow)
    library(broom)
    library(cmprsk)
    library(data.table)
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
out_dir <- file.path("output", "contemporary_6yr_waitlist_pollution_finegray_crr_cengroup_point")
model_result_dir <- file.path(out_dir, "model_results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(model_result_dir, recursive = TRUE, showWarnings = FALSE)

target_organs <- c("HR", "KI", "LI", "LU")
organ_labels <- c(HR = "Heart", KI = "Kidney", LI = "Liver", LU = "Lung")

pollutant_specs <- list(
  pm25 = list(
    years = 2016:2021,
    resolution = "daily",
    dir = file.path("data", "release", "lghap_pm25_zcta_daily_parquet"),
    file = function(year) file.path("data", "release", "lghap_pm25_zcta_daily_parquet", sprintf("lghap_pm25_zcta_daily_%04d.parquet", year)),
    value_col = "pm25_ug_m3",
    term = "pm25_waitlist_5ug",
    raw_col = "pm25_waitlist_ug_m3",
    label = "Daily PM2.5 waitlist-period mean per 5 ug/m3",
    scale = 5
  ),
  o3 = list(
    years = 2019:2024,
    resolution = "daily",
    dir = file.path("data", "release", "o3_zcta_daily_parquet"),
    file = function(year) file.path("data", "release", "o3_zcta_daily_parquet", sprintf("o3_zcta_daily_%04d.parquet", year)),
    value_col = "o3_ppb",
    term = "o3_waitlist_10ppb",
    raw_col = "o3_waitlist_ppb",
    label = "Daily ozone waitlist-period mean per 10 ppb",
    scale = 10
  ),
  no2 = list(
    years = 2020:2025,
    resolution = "monthly",
    dir = file.path("data", "release", "no2_zcta_monthly_parquet"),
    file = function(year) file.path("data", "release", "no2_zcta_monthly_parquet", sprintf("no2_zcta_monthly_%04d.parquet", year)),
    value_col = "no2_ppbv",
    term = "no2_waitlist_10ppb",
    raw_col = "no2_waitlist_ppbv",
    label = "Monthly NO2 waitlist-period mean per 10 ppbv",
    scale = 10
  )
)

progress_path <- file.path(out_dir, "finegray_progress.csv")

log_msg <- function(...) {
  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
  flush.console()
}

append_progress <- function(row) {
  write_csv(as_tibble(row), progress_path, append = file.exists(progress_path))
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
  attr(pb, "total") <- total
  pb
}

tick_progress <- function(pb, value, detail = "") {
  setTxtProgressBar(pb, value)
  if (nzchar(detail)) message("\n", attr(pb, "label"), " ", value, "/", attr(pb, "total"), " | ", detail)
  invisible(NULL)
}

prepare_pollutant_cohort <- function(base, spec) {
  start_date <- as.Date(sprintf("%04d-01-01", min(spec$years)))
  end_date <- as.Date(sprintf("%04d-12-31", max(spec$years)))
  base %>%
    filter(index_date <= end_date, index_date >= start_date) %>%
    mutate(
      pollutant_followup_end_date = pmin(observed_end_date, end_date),
      pollutant_followup_days = as.numeric(pollutant_followup_end_date - index_date),
      fstatus = case_when(
        observed_end_date <= end_date & event_type == "death_or_deterioration_delist" ~ "adverse",
        observed_end_date <= end_date & event_type == "transplant_or_improvement" ~ "transplant_or_improvement",
        TRUE ~ "censor"
      )
    ) %>%
    filter(pollutant_followup_days > 0)
}

compute_daily_exposure <- function(cohort, spec, pollutant_name) {
  cohort_dt <- as.data.table(cohort)
  parts <- list()
  value_col <- spec$value_col

  for (yr in spec$years) {
    log_msg("Computing daily ", pollutant_name, " exposure for ", yr)
    year_start <- as.Date(sprintf("%04d-01-01", yr))
    year_end <- as.Date(sprintf("%04d-12-31", yr))
    intervals_y <- cohort_dt[index_date <= year_end & pollutant_followup_end_date >= year_start]
    if (nrow(intervals_y) == 0L) next
    intervals_y[, `:=`(
      interval_start_date = pmax(index_date, year_start),
      interval_end_date = pmin(pollutant_followup_end_date, year_end)
    )]
    intervals_y <- intervals_y[interval_end_date >= interval_start_date]
    if (nrow(intervals_y) == 0L) next

    vals <- read_parquet(spec$file(yr), col_select = c("zip", "date", value_col)) %>%
      filter(zip %in% unique(intervals_y$candidate_zip)) %>%
      as.data.table()
    setnames(vals, value_col, "pollution_value")
    setorder(vals, zip, date)
    vals[, `:=`(
      valid_value = as.integer(!is.na(pollution_value)),
      value_for_sum = fifelse(is.na(pollution_value), 0, pollution_value)
    )]
    vals[, `:=`(
      cum_value = cumsum(value_for_sum),
      cum_n = cumsum(valid_value)
    ), by = zip]
    cum_lookup <- vals[, .(zip, date, cum_value, cum_n)]

    end_lookup <- intervals_y[, .(row_id = .I, waitlist_row_id, zip = candidate_zip, date = interval_end_date)]
    start_lookup <- intervals_y[, .(row_id = .I, zip = candidate_zip, date = interval_start_date - 1L)]
    end_cum <- cum_lookup[end_lookup, on = c("zip", "date")]
    start_cum <- cum_lookup[start_lookup, on = c("zip", "date")]
    for (col in c("cum_value", "cum_n")) {
      set(end_cum, which(is.na(end_cum[[col]])), col, 0)
      set(start_cum, which(is.na(start_cum[[col]])), col, 0)
    }
    parts[[as.character(yr)]] <- data.table(
      waitlist_row_id = intervals_y$waitlist_row_id,
      exposure_sum = end_cum$cum_value - start_cum$cum_value,
      exposure_days = end_cum$cum_n - start_cum$cum_n
    )
    rm(vals, cum_lookup, end_lookup, start_lookup, end_cum, start_cum)
    invisible(gc())
  }

  rbindlist(parts, use.names = TRUE, fill = TRUE)[
    ,
    .(
      exposure_days = sum(exposure_days, na.rm = TRUE),
      exposure_mean = sum(exposure_sum, na.rm = TRUE) / sum(exposure_days, na.rm = TRUE)
    ),
    by = waitlist_row_id
  ]
}

make_month_intervals <- function(cohort) {
  intervals <- as.data.table(cohort)[, .(waitlist_row_id, candidate_zip, index_date, pollutant_followup_end_date)]
  intervals[, month_seq := lapply(seq_len(.N), function(i) {
    start_month <- as.Date(format(index_date[i], "%Y-%m-01"))
    end_month <- as.Date(format(pollutant_followup_end_date[i], "%Y-%m-01"))
    if (is.na(start_month) || is.na(end_month) || end_month < start_month) return(as.Date(character()))
    seq(start_month, end_month, by = "month")
  })]
  intervals <- intervals[lengths(month_seq) > 0L]
  interval_lengths <- lengths(intervals$month_seq)
  out <- intervals[rep(seq_len(nrow(intervals)), interval_lengths)]
  out[, month_start := as.Date(unlist(intervals$month_seq, use.names = FALSE), origin = "1970-01-01")]
  out[, month_seq := NULL]
  out[, month_end := as.Date(format(month_start + 35L, "%Y-%m-01")) - 1L]
  out[, interval_start_date := pmax(index_date, month_start)]
  out[, interval_end_date := pmin(pollutant_followup_end_date, month_end)]
  out <- out[interval_end_date >= interval_start_date]
  out[, `:=`(
    year = as.integer(format(month_start, "%Y")),
    month = as.integer(format(month_start, "%m")),
    interval_days = as.integer(interval_end_date - interval_start_date) + 1L
  )]
  out
}

compute_monthly_exposure <- function(cohort, spec, pollutant_name) {
  log_msg("Computing monthly ", pollutant_name, " exposure")
  month_intervals <- make_month_intervals(cohort)
  values <- rbindlist(lapply(spec$years, function(yr) {
    read_parquet(spec$file(yr), col_select = c("zip", "year", "month", spec$value_col)) %>%
      as.data.table()
  }), use.names = TRUE)
  setnames(values, spec$value_col, "pollution_value")
  joined <- values[month_intervals, on = c("zip" = "candidate_zip", "year", "month")]
  joined[!is.na(pollution_value), .(
    exposure_days = sum(interval_days),
    exposure_mean = sum(pollution_value * interval_days) / sum(interval_days)
  ), by = waitlist_row_id]
}

fit_finegray <- function(dat, pollutant_name, spec, org) {
  result_path <- file.path(model_result_dir, paste0(tolower(org), "_", pollutant_name, "_finegray.csv"))
  if (file.exists(result_path) && Sys.getenv("FG_OVERWRITE", "FALSE") != "TRUE") {
    log_msg("Skipping existing Fine-Gray result ", org, " | ", pollutant_name)
    return(result_path)
  }

  vars_needed <- c(spec$term, "pollutant_followup_days", "fstatus", "age", "sex", "race", "listing_year", "organ_score", "listing_center")
  model_dat <- dat %>%
    filter(WL_ORG == org, complete.cases(across(all_of(vars_needed))), pollutant_followup_days > 0) %>%
    mutate(
      .fg_time = pmax(pollutant_followup_days, 0.5),
      .fg_status = factor(fstatus, levels = c("censor", "adverse", "transplant_or_improvement")),
      .fg_code = case_when(
        fstatus == "adverse" ~ 1L,
        fstatus == "transplant_or_improvement" ~ 2L,
        TRUE ~ 0L
      ),
      .exposure = .data[[spec$term]],
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

  log_msg("Fine-Gray start ", org, " | ", pollutant_name, " n=", source_n, " adverse=", source_adverse, " competing=", source_competing)
  started_at <- Sys.time()
  append_progress(tibble(
    pollutant = pollutant_name,
    organ = org,
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

  fit <- NULL
  fit_error <- NULL
  cov_mat <- NULL
  elapsed <- system.time({
    cov_mat <- tryCatch(
      model.matrix(
        ~ .exposure + age + sex + race + listing_year + organ_score,
        data = model_dat
      )[, -1, drop = FALSE],
      error = function(e) {
        fit_error <<- conditionMessage(e)
        NULL
      }
    )
    if (!is.null(cov_mat)) {
      keep_cols <- apply(cov_mat, 2, function(x) length(unique(x)) > 1L)
      cov_mat <- cov_mat[, keep_cols, drop = FALSE]
      if (!".exposure" %in% colnames(cov_mat)) {
        fit_error <<- "Exposure column was dropped from the model matrix."
      } else {
        fit <- tryCatch(
          crr(
            ftime = model_dat$.fg_time,
            fstatus = model_dat$.fg_code,
            cov1 = cov_mat,
            failcode = 1,
            cencode = 0,
            cengroup = model_dat$listing_center,
            variance = FALSE,
            maxiter = 100
          ),
          error = function(e) {
            fit_error <<- conditionMessage(e)
            NULL
          }
        )
      }
    }
  })
  finished_at <- Sys.time()

  if (is.null(fit)) {
    result <- tibble(
      organ = org,
      organ_label = recode(org, !!!organ_labels),
      pollutant = pollutant_name,
      exposure_label = spec$label,
      exposure_years = paste(range(spec$years), collapse = "-"),
      endpoint = "death_or_deterioration_delist_subdistribution",
      n = source_n,
      adverse_events = source_adverse,
      competing_transplant_or_improvement = source_competing,
      expanded_rows = NA_real_,
      subdistribution_hazard_ratio = NA_real_,
      conf_low = NA_real_,
      conf_high = NA_real_,
      p_value = NA_real_,
      elapsed_seconds = unname(elapsed[["elapsed"]]),
      adjustment_set = "age + sex + race + listing_year + baseline_organ_score; listing_center used as cengroup for censoring distribution",
      variance_estimator = "cmprsk_crr_model_based_center_cengroup",
      sampled = FALSE,
      error = safe_error(fit_error)
    )
  } else {
    exposure_idx <- which(names(fit$coef) == ".exposure")
    exposure_coef <- fit$coef[[exposure_idx]]
    exposure_se <- NA_real_
    exposure_z <- NA_real_
    result <- tibble(
      organ = org,
      organ_label = recode(org, !!!organ_labels),
      pollutant = pollutant_name,
      exposure_label = spec$label,
      exposure_years = paste(range(spec$years), collapse = "-"),
      endpoint = "death_or_deterioration_delist_subdistribution",
      n = source_n,
      adverse_events = source_adverse,
      competing_transplant_or_improvement = source_competing,
      expanded_rows = NA_real_,
      subdistribution_hazard_ratio = exp(exposure_coef),
      conf_low = NA_real_,
      conf_high = NA_real_,
      p_value = NA_real_,
      elapsed_seconds = unname(elapsed[["elapsed"]]),
      adjustment_set = "age + sex + race + listing_year + baseline_organ_score; listing_center used as cengroup for censoring distribution",
      variance_estimator = "not_estimated_cmprsk_crr_variance_false_center_cengroup",
      sampled = FALSE,
      error = NA_character_
    )
  }
  write_csv(result, result_path)
  append_progress(tibble(
    pollutant = pollutant_name,
    organ = org,
    status = if (is.null(fit)) "error" else "completed",
    started_at = as.character(started_at),
    finished_at = as.character(finished_at),
    elapsed_seconds = unname(elapsed[["elapsed"]]),
    n = source_n,
    adverse_events = source_adverse,
    competing_events = source_competing,
    expanded_rows = NA_real_,
    error = safe_error(fit_error)
  ))
  log_msg("Fine-Gray done ", org, " | ", pollutant_name, " elapsed_seconds=", round(unname(elapsed[["elapsed"]]), 1))
  result_path
}

log_msg("Reading waitlist cohort from ", in_path)
base <- read_csv(in_path, show_col_types = FALSE) %>%
  mutate(
    index_date = as.Date(index_date),
    observed_end_date = as.Date(observed_end_date),
    listing_year = as.integer(as.character(listing_year)),
    sex = factor(sex),
    race = factor(race),
    listing_center = factor(listing_center)
  )

pb <- make_progress(length(pollutant_specs) * length(target_organs), "Contemporary Fine-Gray")
result_paths <- character()
counter <- 0L

for (pollutant_name in names(pollutant_specs)) {
  spec <- pollutant_specs[[pollutant_name]]
  log_msg("Preparing ", pollutant_name, " cohort for years ", paste(range(spec$years), collapse = "-"))
  pollutant_cohort <- prepare_pollutant_cohort(base, spec)
  exposure <- if (spec$resolution == "daily") {
    compute_daily_exposure(pollutant_cohort, spec, pollutant_name)
  } else {
    compute_monthly_exposure(pollutant_cohort, spec, pollutant_name)
  }
  exposure_tbl <- as_tibble(exposure)
  days_col <- paste0(pollutant_name, "_days")
  names(exposure_tbl)[names(exposure_tbl) == "exposure_mean"] <- spec$raw_col
  names(exposure_tbl)[names(exposure_tbl) == "exposure_days"] <- days_col
  if (!all(c("waitlist_row_id", spec$raw_col, days_col) %in% names(exposure_tbl))) {
    stop(
      "Exposure summary is missing expected columns for ", pollutant_name,
      ". Found: ", paste(names(exposure_tbl), collapse = ", ")
    )
  }
  analysis_dat <- pollutant_cohort %>%
    select(-any_of(c(spec$raw_col, days_col, spec$term))) %>%
    left_join(exposure_tbl, by = "waitlist_row_id")
  if (!spec$raw_col %in% names(analysis_dat)) {
    stop("Joined analysis data is missing exposure column: ", spec$raw_col)
  }
  analysis_dat[[spec$term]] <- analysis_dat[[spec$raw_col]] / spec$scale

  write_csv(
    analysis_dat %>%
      select(
        waitlist_row_id, PERS_ID, PX_ID, WL_ORG, candidate_zip, index_date,
        observed_end_date, pollutant_followup_end_date, pollutant_followup_days,
        fstatus, event_type, age, sex, race, listing_year, organ_score,
        organ_score_name, listing_center, all_of(spec$raw_col), all_of(spec$term)
      ),
    file.path(out_dir, paste0(pollutant_name, "_analysis_dataset.csv.gz"))
  )

  for (org in target_organs) {
    counter <- counter + 1L
    tick_progress(pb, counter - 1L, paste(pollutant_name, org, "starting"))
    result_paths <- c(result_paths, fit_finegray(analysis_dat, pollutant_name, spec, org))
    tick_progress(pb, counter, paste(pollutant_name, org, "finished"))
    invisible(gc())
  }
}
close(pb)

results <- bind_rows(lapply(result_paths[file.exists(result_paths)], read_csv, show_col_types = FALSE))
write_csv(results, file.path(out_dir, "contemporary_6yr_waitlist_pollution_finegray_results.csv"))
log_msg("Wrote contemporary six-year Fine-Gray outputs to ", out_dir)
