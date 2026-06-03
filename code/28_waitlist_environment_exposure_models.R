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
  library(data.table)
  library(dplyr)
  library(haven)
  library(readr)
  library(stringr)
  library(survival)
  library(tibble)
})

saf_paths_source <- file.path("code", "saf_paths.R")
if (!file.exists(saf_paths_source)) saf_paths_source <- "saf_paths.R"
source(saf_paths_source)
saf_paths <- get_saf_paths()
assert_saf_files(saf_paths, include_stathist = TRUE)

pubsaf_dir <- saf_paths$pubsaf_dir
supp_dir <- saf_paths$supp_dir
gridmet_dir <- file.path("data", "release", "gridmet_zcta_daily_parquet")
pollution_release_dir <- file.path("data", "release", "air_pollution_zcta_parquet")
out_dir <- file.path("output", "waitlist_environment_exposure")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

analysis_start_year <- as.integer(Sys.getenv("ENV_ANALYSIS_START_YEAR", "2006"))
analysis_end_year <- as.integer(Sys.getenv("ENV_ANALYSIS_END_YEAR", "2023"))
analysis_end_date <- as.Date(sprintf("%04d-12-31", analysis_end_year))
target_organs <- c("HR", "KI", "LI", "LU")

log_msg <- function(...) {
  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
  flush.console()
}

make_progress <- function(total, label) {
  if (total <= 0L) return(NULL)
  pb <- txtProgressBar(min = 0, max = total, style = 3)
  attr(pb, "label") <- label
  attr(pb, "total") <- total
  pb
}

tick_progress <- function(pb, value, detail = "") {
  if (is.null(pb)) return(invisible(NULL))
  setTxtProgressBar(pb, value)
  if (nzchar(detail)) message("\n", attr(pb, "label"), " ", value, "/", attr(pb, "total"), " | ", detail)
  invisible(NULL)
}

close_progress <- function(pb) {
  if (!is.null(pb)) close(pb)
}

clean_zip <- function(x) {
  z <- str_extract(as.character(x), "[0-9]+")
  z <- if_else(nchar(z) > 5L, substr(z, 1, 5), z)
  if_else(is.na(z), NA_character_, str_pad(z, 5, side = "left", pad = "0"))
}

as_factor_missing <- function(x) {
  factor(if_else(is.na(x) | as.character(x) == "", "Missing", as.character(x)))
}

flag_yes <- function(x) {
  y <- str_to_upper(str_trim(as.character(x)))
  as.integer(y %in% c("1", "Y", "YES", "TRUE", "T"))
}

safe_log1p <- function(x) {
  log(pmax(as.numeric(x), 0) + 1)
}

median_impute <- function(x) {
  x <- as.numeric(x)
  med <- median(x, na.rm = TRUE)
  if (!is.finite(med)) med <- NA_real_
  if_else(is.na(x), med, x)
}

collapse_rare <- function(x, min_n = 100L) {
  y <- as.character(x)
  y[is.na(y) | y == ""] <- "Missing"
  tab <- table(y)
  y[!y %in% names(tab)[tab >= min_n]] <- "Other"
  factor(y)
}

calc_egfr_2021 <- function(creatinine, age, sex) {
  scr <- as.numeric(creatinine)
  age <- as.numeric(age)
  sex_chr <- str_to_upper(as.character(sex))
  female <- sex_chr %in% c("F", "FEMALE")
  kappa <- if_else(female, 0.7, 0.9)
  alpha <- if_else(female, -0.241, -0.302)
  sex_mult <- if_else(female, 1.012, 1.0)
  142 * pmin(scr / kappa, 1)^alpha * pmax(scr / kappa, 1)^(-1.200) * 0.9938^age * sex_mult
}

extract_formula_vars <- function(terms) {
  vars <- unique(unlist(str_extract_all(terms, "[A-Za-z][A-Za-z0-9_]*")))
  setdiff(vars, c("strata", "cluster", "frailty", "distribution"))
}

read_pollution_tables <- function() {
  log_msg("Reading air-pollution Parquet release assets")
  pm25 <- read_parquet(file.path(pollution_release_dir, "air_pollution_zcta_pm25_monthly_2005_2023.parquet")) %>%
    transmute(zip = clean_zip(zip), year, month, pm25_ug_m3)
  o3 <- read_parquet(file.path(pollution_release_dir, "air_pollution_zcta_o3_monthly_2005_2023.parquet")) %>%
    transmute(zip = clean_zip(zip), year, month, o3_ppb)
  no2 <- read_parquet(file.path(pollution_release_dir, "air_pollution_zcta_no2_annual_2005_2024.parquet")) %>%
    transmute(zip = clean_zip(zip), year, no2)
  list(pm25 = as.data.table(pm25), o3 = as.data.table(o3), no2 = as.data.table(no2))
}

candidate_files <- tribble(
  ~candidate_group, ~path,
  "kidney_pancreas", file.path(pubsaf_dir, "cand_kipa.sas7bdat"),
  "liver_intestine", file.path(pubsaf_dir, "cand_liin.sas7bdat"),
  "thoracic", file.path(pubsaf_dir, "cand_thor.sas7bdat")
)

candidate_cols <- c(
  "PERS_ID", "PX_ID", "WL_ORG", "CAN_LISTING_DT", "CAN_ACTIVATE_DT", "CAN_SOURCE",
  "CAN_REM_DT", "CAN_REM_CD", "CAN_DEATH_DT", "CAN_ENDWLFU", "CAN_AGE_AT_LISTING",
  "CAN_GENDER", "CAN_RACE_SRTR", "CAN_LISTING_CTR_CD", "CAN_INIT_STAT",
  "CAN_ON_DIAL", "CAN_DIAL", "CAN_DIAL_DT", "CAN_MOST_RECENT_CREAT",
  "CAN_INIT_SRTR_LAB_MELD", "CAN_LAST_SRTR_LAB_MELD", "CAN_TOT_BILI",
  "CAN_TOT_ALBUMIN", "CAN_LAST_SERUM_SODIUM", "CAN_FVC", "CAN_FEV1", "CAN_PCO2",
  "CAN_AT_REST_O2", "CAN_SIX_MIN_WALK", "CAN_PULM_ART_MEAN", "CAN_CARDIAC_OUTPUT",
  "CAN_VENTILATOR", "CAN_ON_VENTILATOR", "CAN_ECMO", "CAN_VAD_TAH", "CAN_VAD1",
  "CAN_VAD2", "CAN_CORTICOST_DEPND"
)

log_msg("Reading SAF candidate and ZIP inputs")
candidate_zip <- read_sas(file.path(supp_dir, "canzip2506.sas7bdat")) %>%
  transmute(PERS_ID, PX_ID, candidate_zip = clean_zip(CAN_PERM_ZIP)) %>%
  distinct(PERS_ID, PX_ID, .keep_all = TRUE)

candidate_data <- candidate_files %>%
  rowwise() %>%
  mutate(data = list(read_sas(path, col_select = any_of(candidate_cols)) %>% mutate(candidate_group = candidate_group))) %>%
  ungroup() %>%
  select(data) %>%
  tidyr::unnest(data)

heart_just <- read_sas(
  file.path(pubsaf_dir, "statjust_hr1a.sas7bdat"),
  col_select = any_of(c(
    "PX_ID", "WL_ORG", "CANHX_CHG_DT", "CANHX_LVAD_TYPE", "CANHX_RVAD_TYPE",
    "CANHX_ECMO", "CANHX_LAB_SERUM_CREAT", "CANHX_LAB_BILI",
    "CANHX_LAB_ALBUMIN", "CANHX_LAB_SODIUM", "CANHX_LAB_BNP"
  ))
) %>%
  filter(WL_ORG == "HR") %>%
  arrange(PX_ID, CANHX_CHG_DT) %>%
  group_by(PX_ID) %>%
  slice_tail(n = 1) %>%
  ungroup() %>%
  transmute(
    PX_ID,
    hr_short_mcs = as.integer(flag_yes(CANHX_ECMO) == 1L | flag_yes(CANHX_RVAD_TYPE) == 1L),
    hr_durable_lvad = flag_yes(CANHX_LVAD_TYPE),
    hr_creatinine = CANHX_LAB_SERUM_CREAT,
    hr_bilirubin = CANHX_LAB_BILI,
    hr_albumin = CANHX_LAB_ALBUMIN,
    hr_sodium = CANHX_LAB_SODIUM,
    hr_bnp = CANHX_LAB_BNP
  )

log_msg("Constructing waitlist cohort")
cohort <- candidate_data %>%
  filter(WL_ORG %in% target_organs) %>%
  left_join(candidate_zip, by = c("PERS_ID", "PX_ID")) %>%
  left_join(heart_just, by = "PX_ID") %>%
  mutate(
    waitlist_row_id = row_number(),
    index_date = coalesce(CAN_ACTIVATE_DT, CAN_LISTING_DT),
    index_year = as.integer(format(index_date, "%Y")),
    event_date = coalesce(CAN_REM_DT, CAN_DEATH_DT, CAN_ENDWLFU),
    observed_end_date = pmin(event_date, analysis_end_date),
    event_observed = !is.na(event_date) & event_date <= analysis_end_date,
    adverse_event = as.integer(event_observed & CAN_REM_CD %in% c(8, 13)),
    favorable_exit = as.integer(event_observed & CAN_REM_CD %in% c(4, 12, 14, 15, 18, 19, 21, 22, 23)),
    transplant_event = as.integer(event_observed & CAN_REM_CD %in% c(4, 14, 15, 18, 19, 21, 22, 23)),
    improved_delist = as.integer(event_observed & CAN_REM_CD == 12),
    other_exit = as.integer(event_observed & adverse_event == 0L & favorable_exit == 0L & !is.na(CAN_REM_CD)),
    followup_days = as.numeric(observed_end_date - index_date),
    age = CAN_AGE_AT_LISTING,
    sex = as_factor_missing(if_else(CAN_GENDER %in% c("M", "F"), CAN_GENDER, NA_character_)),
    race = collapse_rare(CAN_RACE_SRTR, min_n = 100L),
    listing_center_raw = as.character(CAN_LISTING_CTR_CD),
    listing_center = collapse_rare(CAN_LISTING_CTR_CD, min_n = 100L),
    index_year_centered = index_year - 2015,
    kidney_on_dialysis = flag_yes(coalesce(as.character(CAN_ON_DIAL), as.character(CAN_DIAL))),
    dialysis_years = as.numeric(index_date - CAN_DIAL_DT) / 365.25,
    dialysis_years = if_else(is.na(dialysis_years) | dialysis_years < 0, 0, dialysis_years),
    kidney_dialysis_age_score = age / 10 + 0.8 * kidney_on_dialysis + 0.3 * log1p(dialysis_years),
    liver_meld = coalesce(CAN_INIT_SRTR_LAB_MELD, CAN_LAST_SRTR_LAB_MELD),
    hr_creatinine_for_score = coalesce(hr_creatinine, CAN_MOST_RECENT_CREAT),
    hr_albumin = coalesce(hr_albumin, CAN_TOT_ALBUMIN),
    hr_bilirubin = coalesce(hr_bilirubin, CAN_TOT_BILI),
    hr_sodium = coalesce(hr_sodium, CAN_LAST_SERUM_SODIUM),
    hr_short_mcs = coalesce(hr_short_mcs, flag_yes(CAN_ECMO), 0L),
    hr_durable_lvad = coalesce(hr_durable_lvad, flag_yes(CAN_VAD_TAH), 0L)
  ) %>%
  group_by(WL_ORG) %>%
  mutate(
    hr_creatinine_imputed = if_else(WL_ORG == "HR", median_impute(hr_creatinine_for_score), as.numeric(hr_creatinine_for_score)),
    hr_bilirubin_imputed = if_else(WL_ORG == "HR", median_impute(hr_bilirubin), as.numeric(hr_bilirubin)),
    hr_albumin_imputed = if_else(WL_ORG == "HR", median_impute(hr_albumin), as.numeric(hr_albumin)),
    hr_sodium_imputed = if_else(WL_ORG == "HR", median_impute(hr_sodium), as.numeric(hr_sodium)),
    hr_bnp_imputed = if_else(WL_ORG == "HR", median_impute(hr_bnp), as.numeric(hr_bnp))
  ) %>%
  ungroup() %>%
  mutate(
    hr_egfr = calc_egfr_2021(hr_creatinine_imputed, age, sex),
    us_crs_proxy = -0.656 * hr_albumin_imputed +
      0.617 * safe_log1p(hr_bilirubin_imputed) -
      0.012 * hr_egfr -
      0.077 * hr_sodium_imputed -
      0.377 * hr_durable_lvad +
      1.092 * hr_short_mcs +
      0.433 * safe_log1p(hr_bnp_imputed),
    lung_las_cas_component_score = rowSums(cbind(
      as.numeric(coalesce(CAN_FEV1, 0)) * -0.01,
      as.numeric(coalesce(CAN_FVC, 0)) * -0.005,
      as.numeric(coalesce(CAN_PCO2, 0)) * 0.02,
      as.numeric(coalesce(CAN_AT_REST_O2, 0)) * 0.01,
      as.numeric(coalesce(CAN_SIX_MIN_WALK, 0)) * -0.002,
      as.numeric(coalesce(CAN_PULM_ART_MEAN, 0)) * 0.01,
      as.numeric(coalesce(CAN_CARDIAC_OUTPUT, 0)) * -0.05,
      flag_yes(coalesce(as.character(CAN_VENTILATOR), as.character(CAN_ON_VENTILATOR))) * 1.0,
      flag_yes(CAN_ECMO) * 1.5,
      flag_yes(CAN_CORTICOST_DEPND) * 0.3
    ), na.rm = TRUE),
    organ_score = case_when(
      WL_ORG == "HR" ~ us_crs_proxy,
      WL_ORG == "KI" ~ kidney_dialysis_age_score,
      WL_ORG == "LI" ~ liver_meld,
      WL_ORG == "LU" ~ lung_las_cas_component_score
    )
  ) %>%
  filter(
    index_year >= analysis_start_year,
    index_year <= analysis_end_year,
    CAN_SOURCE %in% c("A", "R"),
    !is.na(index_date),
    !is.na(observed_end_date),
    !is.na(followup_days),
    followup_days >= 0,
    observed_end_date >= index_date,
    !is.na(candidate_zip)
  ) %>%
  mutate(followup_days = pmax(followup_days, 0.5))

write_csv(
  cohort %>%
    group_by(WL_ORG) %>%
    summarise(
      n = n(),
      people = n_distinct(PERS_ID),
      adverse_events = sum(adverse_event),
      favorable_exits = sum(favorable_exit),
      transplants = sum(transplant_event),
      improved_delist = sum(improved_delist),
      other_exits = sum(other_exit),
      censored_at_environment_end = sum(!event_observed),
      median_followup_days = median(followup_days),
      .groups = "drop"
    ),
  file.path(out_dir, "waitlist_environment_cohort_summary.csv")
)

make_year_intervals <- function(cohort_dt, years) {
  intervals <- cohort_dt[, .(
    waitlist_row_id, candidate_zip, index_date, observed_end_date
  )]
  intervals[, interval_years := lapply(seq_len(.N), function(i) {
    start_year <- max(as.integer(format(index_date[i], "%Y")), min(years))
    stop_year <- min(as.integer(format(observed_end_date[i], "%Y")), max(years))
    if (is.na(start_year) || is.na(stop_year) || stop_year < start_year) return(integer())
    seq.int(start_year, stop_year)
  })]
  intervals <- intervals[lengths(interval_years) > 0L]
  interval_lengths <- lengths(intervals$interval_years)
  out <- intervals[rep(seq_len(nrow(intervals)), interval_lengths)]
  out[, interval_year := unlist(intervals$interval_years, use.names = FALSE)]
  out[, interval_years := NULL]
  out[, interval_start_date := pmax(index_date, as.Date(paste0(interval_year, "-01-01")))]
  out[, interval_end_date := pmin(observed_end_date, as.Date(paste0(interval_year, "-12-31")))]
  out <- out[interval_end_date >= interval_start_date]
  out[, interval_days := as.integer(interval_end_date - interval_start_date) + 1L]
  out
}

make_month_intervals <- function(cohort_dt) {
  intervals <- cohort_dt[, .(waitlist_row_id, candidate_zip, index_date, observed_end_date)]
  intervals[, month_seq := lapply(seq_len(.N), function(i) {
    start_month <- as.Date(format(index_date[i], "%Y-%m-01"))
    end_month <- as.Date(format(observed_end_date[i], "%Y-%m-01"))
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
  out[, interval_end_date := pmin(observed_end_date, month_end)]
  out <- out[interval_end_date >= interval_start_date]
  out[, `:=`(
    year = as.integer(format(month_start, "%Y")),
    month = as.integer(format(month_start, "%m")),
    interval_days = as.integer(interval_end_date - interval_start_date) + 1L
  )]
  out
}

cohort_dt <- as.data.table(cohort)
year_intervals <- make_year_intervals(cohort_dt, analysis_start_year:analysis_end_year)
month_intervals <- make_month_intervals(cohort_dt)

gridmet_vars <- c("tmax_c", "tmin_c", "tmean_c", "temp_range_c", "rmax_pct", "rmin_pct", "rhmean_pct", "sph")
gridmet_cache_dir <- file.path(out_dir, "gridmet_zcta_monthly_cache")
dir.create(gridmet_cache_dir, recursive = TRUE, showWarnings = FALSE)

log_msg("Aggregating daily gridMET to ZCTA-month caches and waitlist means")
gridmet_results <- list()
pb <- make_progress(length(analysis_start_year:analysis_end_year), "gridMET yearly aggregation")
for (i in seq_along(analysis_start_year:analysis_end_year)) {
  yr <- (analysis_start_year:analysis_end_year)[[i]]
  tick_progress(pb, i, as.character(yr))
  intervals_y <- month_intervals[year == yr]
  if (nrow(intervals_y) == 0L) next

  cache_path <- file.path(gridmet_cache_dir, sprintf("gridmet_zcta_monthly_%04d.csv.gz", yr))
  if (file.exists(cache_path)) {
    gm_month <- read_csv(cache_path, col_types = cols(zip = col_character(), .default = col_double()), show_col_types = FALSE) %>%
      as.data.table()
  } else {
    gm_path <- file.path(gridmet_dir, sprintf("gridmet_zcta_daily_%04d.parquet", yr))
    gm_month <- read_parquet(gm_path, col_select = c("zip", "date", gridmet_vars)) %>%
      mutate(
        zip = clean_zip(zip),
        year = as.integer(format(date, "%Y")),
        month = as.integer(format(date, "%m"))
      ) %>%
      filter(zip %in% unique(intervals_y$candidate_zip)) %>%
      group_by(zip, year, month) %>%
      summarise(across(all_of(gridmet_vars), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
    write_csv(gm_month, cache_path)
    gm_month <- as.data.table(gm_month)
  }

  gm_join <- gm_month[
    intervals_y,
    on = c("zip" = "candidate_zip", "year", "month")
  ]
  for (v in gridmet_vars) gm_join[, paste0(v, "_weighted") := get(v) * interval_days]
  gridmet_results[[as.character(yr)]] <- gm_join[
    !is.na(tmean_c),
    c("waitlist_row_id", "interval_days", paste0(gridmet_vars, "_weighted")),
    with = FALSE
  ]
  rm(gm_month, gm_join)
  invisible(gc())
}
close_progress(pb)

gridmet_exposure <- rbindlist(gridmet_results, use.names = TRUE, fill = TRUE)[
  ,
  c(
    list(gridmet_days = sum(interval_days, na.rm = TRUE)),
    lapply(.SD, function(x) sum(x, na.rm = TRUE))
  ),
  by = waitlist_row_id,
  .SDcols = paste0(gridmet_vars, "_weighted")
]
for (v in gridmet_vars) {
  gridmet_exposure[, paste0(v, "_waitlist_mean") := get(paste0(v, "_weighted")) / gridmet_days]
}
gridmet_exposure <- gridmet_exposure[, c("waitlist_row_id", "gridmet_days", paste0(gridmet_vars, "_waitlist_mean")), with = FALSE]

pollution <- read_pollution_tables()

log_msg("Aggregating monthly PM2.5 and ozone over each observed waitlist interval")
pm25_month <- pollution$pm25[month_intervals, on = c("zip" = "candidate_zip", "year", "month")]
o3_month <- pollution$o3[month_intervals, on = c("zip" = "candidate_zip", "year", "month")]
pm25_exp <- pm25_month[!is.na(pm25_ug_m3), .(
  pm25_days = sum(interval_days),
  pm25_waitlist_ug_m3 = sum(pm25_ug_m3 * interval_days) / sum(interval_days)
), by = waitlist_row_id]
o3_exp <- o3_month[!is.na(o3_ppb), .(
  o3_days = sum(interval_days),
  o3_waitlist_ppb = sum(o3_ppb * interval_days) / sum(interval_days)
), by = waitlist_row_id]

log_msg("Aggregating annual NO2 over each observed waitlist interval")
no2_year <- pollution$no2[year_intervals, on = c("zip" = "candidate_zip", "year" = "interval_year")]
no2_exp <- no2_year[!is.na(no2), .(
  no2_days = sum(interval_days),
  no2_waitlist = sum(no2 * interval_days) / sum(interval_days)
), by = waitlist_row_id]

exposure_cohort <- cohort %>%
  left_join(as_tibble(gridmet_exposure), by = "waitlist_row_id") %>%
  left_join(as_tibble(pm25_exp), by = "waitlist_row_id") %>%
  left_join(as_tibble(o3_exp), by = "waitlist_row_id") %>%
  left_join(as_tibble(no2_exp), by = "waitlist_row_id") %>%
  mutate(
    tmean_waitlist_5c = tmean_c_waitlist_mean / 5,
    tmax_waitlist_5c = tmax_c_waitlist_mean / 5,
    tmin_waitlist_5c = tmin_c_waitlist_mean / 5,
    temp_range_waitlist_5c = temp_range_c_waitlist_mean / 5,
    rhmean_waitlist_10pct = rhmean_pct_waitlist_mean / 10,
    rmax_waitlist_10pct = rmax_pct_waitlist_mean / 10,
    rmin_waitlist_10pct = rmin_pct_waitlist_mean / 10,
    pm25_waitlist_5ug = pm25_waitlist_ug_m3 / 5,
    o3_waitlist_10ppb = o3_waitlist_ppb / 10,
    no2_waitlist_10unit = no2_waitlist / 10,
    adverse_vs_favorable = case_when(
      adverse_event == 1L ~ 1L,
      favorable_exit == 1L ~ 0L,
      TRUE ~ NA_integer_
    )
  )

write_csv(
  exposure_cohort %>%
    group_by(WL_ORG) %>%
    summarise(
      n = n(),
      complete_environment = sum(complete.cases(
        tmean_waitlist_5c, temp_range_waitlist_5c, rhmean_waitlist_10pct,
        pm25_waitlist_5ug, o3_waitlist_10ppb, no2_waitlist_10unit
      )),
      adverse_events = sum(adverse_event),
      favorable_exits = sum(favorable_exit),
      adverse_vs_favorable_complete = sum(!is.na(adverse_vs_favorable) & complete.cases(
        tmean_waitlist_5c, temp_range_waitlist_5c, rhmean_waitlist_10pct,
        pm25_waitlist_5ug, o3_waitlist_10ppb, no2_waitlist_10unit
      )),
      mean_tmean_c = mean(tmean_c_waitlist_mean, na.rm = TRUE),
      mean_rhmean_pct = mean(rhmean_pct_waitlist_mean, na.rm = TRUE),
      mean_pm25 = mean(pm25_waitlist_ug_m3, na.rm = TRUE),
      mean_o3 = mean(o3_waitlist_ppb, na.rm = TRUE),
      mean_no2 = mean(no2_waitlist, na.rm = TRUE),
      .groups = "drop"
    ),
  file.path(out_dir, "waitlist_environment_exposure_summary.csv")
)

write_csv(
  exposure_cohort %>%
    select(
      waitlist_row_id, PERS_ID, PX_ID, WL_ORG, candidate_zip, index_date, observed_end_date,
      followup_days, adverse_event, favorable_exit, transplant_event, improved_delist, other_exit,
      tmax_c_waitlist_mean, tmin_c_waitlist_mean, tmean_c_waitlist_mean,
      temp_range_c_waitlist_mean, rmax_pct_waitlist_mean, rmin_pct_waitlist_mean,
      rhmean_pct_waitlist_mean, sph_waitlist_mean, pm25_waitlist_ug_m3,
      o3_waitlist_ppb, no2_waitlist
    ),
  file.path(out_dir, "waitlist_environment_exposures.csv.gz")
)

environment_terms <- c(
  "tmean_waitlist_5c",
  "temp_range_waitlist_5c",
  "rhmean_waitlist_10pct",
  "pm25_waitlist_5ug",
  "o3_waitlist_10ppb",
  "no2_waitlist_10unit"
)

term_labels <- c(
  tmean_waitlist_5c = "Mean temperature per 5 C",
  temp_range_waitlist_5c = "Diurnal temperature range per 5 C",
  rhmean_waitlist_10pct = "Mean relative humidity per 10 pct",
  pm25_waitlist_5ug = "PM2.5 per 5 ug/m3",
  o3_waitlist_10ppb = "Ozone per 10 ppb",
  no2_waitlist_10unit = "NO2 per 10 units"
)

model_specs <- list(
  temperature_humidity = c("tmean_waitlist_5c", "temp_range_waitlist_5c", "rhmean_waitlist_10pct"),
  pollution = c("pm25_waitlist_5ug", "o3_waitlist_10ppb", "no2_waitlist_10unit"),
  combined_environment = environment_terms
)

adjustment_terms <- c(
  "age", "sex", "race", "index_year_centered", "organ_score",
  "strata(listing_center)", "cluster(PERS_ID)"
)

fit_cox <- function(dat, org, spec_name, env_terms) {
  rhs <- c(env_terms, adjustment_terms)
  vars_needed <- intersect(extract_formula_vars(rhs), names(dat))
  model_dat <- dat %>%
    filter(complete.cases(across(all_of(vars_needed))), followup_days > 0) %>%
    droplevels()
  if (nrow(model_dat) == 0L || sum(model_dat$adverse_event) < 100L) return(tibble())
  form <- as.formula(paste("Surv(followup_days, adverse_event) ~", paste(rhs, collapse = " + ")))
  log_msg("Fitting Cox ", org, " | ", spec_name, " n=", nrow(model_dat), " adverse=", sum(model_dat$adverse_event))
  fit <- coxph(form, data = model_dat, ties = "efron", robust = TRUE, x = FALSE, y = FALSE)
  tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(term %in% env_terms) %>%
    transmute(
      organ = org,
      model = spec_name,
      endpoint = "death_or_deterioration_delist_cause_specific",
      term,
      exposure = recode(term, !!!term_labels),
      n = nrow(model_dat),
      adverse_events = sum(model_dat$adverse_event),
      hazard_ratio = estimate,
      conf_low = conf.low,
      conf_high = conf.high,
      p_value = p.value,
      concordance = unname(summary(fit)$concordance[1])
    )
}

fit_logistic <- function(dat, org, spec_name, env_terms) {
  rhs <- c(env_terms, "age", "sex", "race", "index_year_centered", "organ_score")
  vars_needed <- intersect(extract_formula_vars(rhs), names(dat))
  model_dat <- dat %>%
    filter(!is.na(adverse_vs_favorable), complete.cases(across(all_of(vars_needed)))) %>%
    droplevels()
  if (nrow(model_dat) == 0L || sum(model_dat$adverse_vs_favorable) < 100L) return(tibble())
  form <- as.formula(paste("adverse_vs_favorable ~", paste(rhs, collapse = " + ")))
  log_msg("Fitting logistic ", org, " | ", spec_name, " n=", nrow(model_dat), " adverse=", sum(model_dat$adverse_vs_favorable))
  fit <- glm(form, data = model_dat, family = binomial())
  tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(term %in% env_terms) %>%
    transmute(
      organ = org,
      model = spec_name,
      endpoint = "death_or_deterioration_delist_vs_transplant_or_improvement",
      term,
      exposure = recode(term, !!!term_labels),
      n = nrow(model_dat),
      adverse_events = sum(model_dat$adverse_vs_favorable),
      favorable_exits = nrow(model_dat) - sum(model_dat$adverse_vs_favorable),
      odds_ratio = estimate,
      conf_low = conf.low,
      conf_high = conf.high,
      p_value = p.value
    )
}

log_msg("Fitting organ-specific environment models")
cox_results <- bind_rows(lapply(target_organs, function(org) {
  dat <- exposure_cohort %>% filter(WL_ORG == org)
  bind_rows(lapply(names(model_specs), function(spec_name) fit_cox(dat, org, spec_name, model_specs[[spec_name]])))
}))
write_csv(cox_results, file.path(out_dir, "waitlist_environment_adverse_cox_results.csv"))

logistic_results <- bind_rows(lapply(target_organs, function(org) {
  dat <- exposure_cohort %>% filter(WL_ORG == org)
  bind_rows(lapply(names(model_specs), function(spec_name) fit_logistic(dat, org, spec_name, model_specs[[spec_name]])))
}))

write_csv(logistic_results, file.path(out_dir, "waitlist_environment_adverse_vs_favorable_logistic_results.csv"))

log_msg("Wrote waitlist temperature/humidity/pollution outputs to ", out_dir)
