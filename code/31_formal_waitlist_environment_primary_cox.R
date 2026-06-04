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
  library(tidyr)
})

saf_paths_source <- file.path("code", "saf_paths.R")
if (!file.exists(saf_paths_source)) saf_paths_source <- "saf_paths.R"
source(saf_paths_source)
saf_paths <- get_saf_paths()
assert_saf_files(saf_paths, include_stathist = TRUE)

pubsaf_dir <- saf_paths$pubsaf_dir
supp_dir <- saf_paths$supp_dir
gridmet_dir <- file.path("data", "release", "gridmet_zcta_daily_parquet")
pollution_dir <- file.path("data", "release", "air_pollution_zcta_parquet")
out_dir <- file.path("output", "formal_waitlist_environment")
cache_dir <- file.path(out_dir, "cache")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

analysis_start_year <- as.integer(Sys.getenv("ENV_ANALYSIS_START_YEAR", "2006"))
analysis_end_year <- as.integer(Sys.getenv("ENV_ANALYSIS_END_YEAR", "2023"))
analysis_end_date <- as.Date(sprintf("%04d-12-31", analysis_end_year))
target_organs <- c("HR", "KI", "LI", "LU")
organ_labels <- c(HR = "Heart", KI = "Kidney", LI = "Liver", LU = "Lung")

log_msg <- function(...) {
  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
  flush.console()
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

make_year_intervals <- function(cohort_dt, years) {
  intervals <- cohort_dt[, .(waitlist_row_id, candidate_zip, index_date, observed_end_date)]
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

log_msg("Reading SAF candidate, ZIP, and baseline score inputs")
candidate_zip <- read_sas(file.path(supp_dir, "canzip2506.sas7bdat")) %>%
  transmute(PERS_ID, PX_ID, candidate_zip = clean_zip(CAN_PERM_ZIP)) %>%
  distinct(PERS_ID, PX_ID, .keep_all = TRUE)

candidate_data <- candidate_files %>%
  rowwise() %>%
  mutate(data = list(read_sas(path, col_select = any_of(candidate_cols)) %>% mutate(candidate_group = candidate_group))) %>%
  ungroup() %>%
  select(data) %>%
  unnest(data)

heart_baseline <- read_sas(
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
  slice_head(n = 1) %>%
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

cohort <- candidate_data %>%
  filter(WL_ORG %in% target_organs) %>%
  left_join(candidate_zip, by = c("PERS_ID", "PX_ID")) %>%
  left_join(heart_baseline, by = "PX_ID") %>%
  mutate(
    waitlist_row_id = row_number(),
    index_date = coalesce(CAN_ACTIVATE_DT, CAN_LISTING_DT),
    index_year = as.integer(format(index_date, "%Y")),
    event_date = coalesce(CAN_REM_DT, CAN_DEATH_DT, CAN_ENDWLFU),
    observed_end_date = pmin(event_date, analysis_end_date),
    event_observed = !is.na(event_date) & event_date <= analysis_end_date,
    adverse_event = as.integer(event_observed & CAN_REM_CD %in% c(8, 13)),
    transplant_or_improvement = as.integer(event_observed & CAN_REM_CD %in% c(4, 12, 14, 15, 18, 19, 21, 22, 23)),
    other_exit = as.integer(event_observed & adverse_event == 0L & transplant_or_improvement == 0L & !is.na(CAN_REM_CD)),
    event_type = factor(
      case_when(
        adverse_event == 1L ~ "adverse",
        transplant_or_improvement == 1L ~ "transplant_or_improvement",
        other_exit == 1L ~ "other_exit",
        TRUE ~ "censor"
      ),
      levels = c("censor", "adverse", "transplant_or_improvement", "other_exit")
    ),
    followup_days = as.numeric(observed_end_date - index_date),
    age = as.numeric(CAN_AGE_AT_LISTING),
    sex = as_factor_missing(if_else(CAN_GENDER %in% c("M", "F"), CAN_GENDER, NA_character_)),
    race = as_factor_missing(CAN_RACE_SRTR),
    listing_center = as_factor_missing(CAN_LISTING_CTR_CD),
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
    ),
    organ_score_name = case_when(
      WL_ORG == "HR" ~ "US-CRS proxy baseline",
      WL_ORG == "KI" ~ "Dialysis + age baseline",
      WL_ORG == "LI" ~ "Initial MELD, last if initial missing",
      WL_ORG == "LU" ~ "LAS/CAS component proxy baseline"
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

cohort_dt <- as.data.table(cohort)
year_intervals <- make_year_intervals(cohort_dt, analysis_start_year:analysis_end_year)
month_intervals <- make_month_intervals(cohort_dt)

log_msg("Computing exact daily tmax/rmax waitlist means from gridMET")
daily_cache_path <- file.path(cache_dir, "waitlist_daily_tmax_rmax_exact.csv.gz")
if (file.exists(daily_cache_path)) {
  daily_exposure <- read_csv(daily_cache_path, show_col_types = FALSE) %>% as.data.table()
} else {
  daily_parts <- list()
  for (yr in analysis_start_year:analysis_end_year) {
    log_msg("Daily gridMET cumulative lookup for ", yr)
    intervals_y <- year_intervals[interval_year == yr]
    if (nrow(intervals_y) == 0L) next
    gm_path <- file.path(gridmet_dir, sprintf("gridmet_zcta_daily_%04d.parquet", yr))
    gm <- read_parquet(gm_path, col_select = c("zip", "date", "tmax_c", "rmax_pct")) %>%
      mutate(zip = clean_zip(zip)) %>%
      filter(zip %in% unique(intervals_y$candidate_zip)) %>%
      as.data.table()
    setorder(gm, zip, date)
    gm[, `:=`(
      valid_tmax = as.integer(!is.na(tmax_c)),
      valid_rmax = as.integer(!is.na(rmax_pct)),
      tmax_for_sum = fifelse(is.na(tmax_c), 0, tmax_c),
      rmax_for_sum = fifelse(is.na(rmax_pct), 0, rmax_pct)
    )]
    gm[, `:=`(
      cum_tmax = cumsum(tmax_for_sum),
      cum_rmax = cumsum(rmax_for_sum),
      cum_n_tmax = cumsum(valid_tmax),
      cum_n_rmax = cumsum(valid_rmax)
    ), by = zip]
    cum_lookup <- gm[, .(zip, date, cum_tmax, cum_rmax, cum_n_tmax, cum_n_rmax)]

    end_lookup <- intervals_y[, .(row_id = .I, waitlist_row_id, zip = candidate_zip, date = interval_end_date)]
    start_lookup <- intervals_y[, .(row_id = .I, zip = candidate_zip, date = interval_start_date - 1L)]
    end_cum <- cum_lookup[end_lookup, on = c("zip", "date")]
    start_cum <- cum_lookup[start_lookup, on = c("zip", "date")]
    for (col in c("cum_tmax", "cum_rmax", "cum_n_tmax", "cum_n_rmax")) {
      set(end_cum, which(is.na(end_cum[[col]])), col, 0)
      set(start_cum, which(is.na(start_cum[[col]])), col, 0)
    }
    daily_parts[[as.character(yr)]] <- data.table(
      waitlist_row_id = intervals_y$waitlist_row_id,
      tmax_sum = end_cum$cum_tmax - start_cum$cum_tmax,
      rmax_sum = end_cum$cum_rmax - start_cum$cum_rmax,
      tmax_days = end_cum$cum_n_tmax - start_cum$cum_n_tmax,
      rmax_days = end_cum$cum_n_rmax - start_cum$cum_n_rmax
    )
    rm(gm, cum_lookup, end_lookup, start_lookup, end_cum, start_cum)
    invisible(gc())
  }
  daily_exposure <- rbindlist(daily_parts, use.names = TRUE, fill = TRUE)[
    ,
    .(
      tmax_days = sum(tmax_days, na.rm = TRUE),
      rmax_days = sum(rmax_days, na.rm = TRUE),
      tmax_c_waitlist_mean_exact_daily = sum(tmax_sum, na.rm = TRUE) / sum(tmax_days, na.rm = TRUE),
      rmax_pct_waitlist_mean_exact_daily = sum(rmax_sum, na.rm = TRUE) / sum(rmax_days, na.rm = TRUE)
    ),
    by = waitlist_row_id
  ]
  write_csv(as_tibble(daily_exposure), daily_cache_path)
}

log_msg("Aggregating monthly/annual air pollution over waitlist intervals")
pm25 <- read_parquet(file.path(pollution_dir, "air_pollution_zcta_pm25_monthly_2005_2023.parquet")) %>%
  transmute(zip = clean_zip(zip), year, month, pm25_ug_m3) %>%
  as.data.table()
o3 <- read_parquet(file.path(pollution_dir, "air_pollution_zcta_o3_monthly_2005_2023.parquet")) %>%
  transmute(zip = clean_zip(zip), year, month, o3_ppb) %>%
  as.data.table()
no2 <- read_parquet(file.path(pollution_dir, "air_pollution_zcta_no2_annual_2005_2024.parquet")) %>%
  transmute(zip = clean_zip(zip), year, no2) %>%
  as.data.table()

pm25_month <- pm25[month_intervals, on = c("zip" = "candidate_zip", "year", "month")]
o3_month <- o3[month_intervals, on = c("zip" = "candidate_zip", "year", "month")]
no2_year <- no2[year_intervals, on = c("zip" = "candidate_zip", "year" = "interval_year")]

pm25_exp <- pm25_month[!is.na(pm25_ug_m3), .(
  pm25_days = sum(interval_days),
  pm25_waitlist_ug_m3 = sum(pm25_ug_m3 * interval_days) / sum(interval_days)
), by = waitlist_row_id]
o3_exp <- o3_month[!is.na(o3_ppb), .(
  o3_days = sum(interval_days),
  o3_waitlist_ppb = sum(o3_ppb * interval_days) / sum(interval_days)
), by = waitlist_row_id]
no2_exp <- no2_year[!is.na(no2), .(
  no2_days = sum(interval_days),
  no2_waitlist = sum(no2 * interval_days) / sum(interval_days)
), by = waitlist_row_id]

analysis_dat <- cohort %>%
  left_join(as_tibble(daily_exposure), by = "waitlist_row_id") %>%
  left_join(as_tibble(pm25_exp), by = "waitlist_row_id") %>%
  left_join(as_tibble(o3_exp), by = "waitlist_row_id") %>%
  left_join(as_tibble(no2_exp), by = "waitlist_row_id") %>%
  mutate(
    tmax_waitlist_5c = tmax_c_waitlist_mean_exact_daily / 5,
    rmax_waitlist_10pct = rmax_pct_waitlist_mean_exact_daily / 10,
    pm25_waitlist_5ug = pm25_waitlist_ug_m3 / 5,
    o3_waitlist_10ppb = o3_waitlist_ppb / 10,
    no2_waitlist_10unit = no2_waitlist / 10
  )

write_csv(
  analysis_dat %>%
    select(
      waitlist_row_id, PERS_ID, PX_ID, WL_ORG, candidate_zip, index_date, observed_end_date,
      followup_days, event_type, adverse_event, transplant_or_improvement, other_exit,
      age, sex, race, index_year, index_year_centered, organ_score, organ_score_name, listing_center,
      tmax_c_waitlist_mean_exact_daily, rmax_pct_waitlist_mean_exact_daily,
      pm25_waitlist_ug_m3, o3_waitlist_ppb, no2_waitlist
    ),
  file.path(out_dir, "formal_waitlist_environment_analysis_dataset.csv.gz")
)

exposure_specs <- list(
  heat_humidity = c("tmax_waitlist_5c", "rmax_waitlist_10pct"),
  pm25 = "pm25_waitlist_5ug",
  o3 = "o3_waitlist_10ppb",
  no2 = "no2_waitlist_10unit"
)

term_labels <- c(
  tmax_waitlist_5c = "Tmax per 5 C",
  rmax_waitlist_10pct = "Maximum relative humidity per 10 pct",
  pm25_waitlist_5ug = "PM2.5 per 5 ug/m3",
  o3_waitlist_10ppb = "Ozone per 10 ppb",
  no2_waitlist_10unit = "NO2 per 10 units"
)

fit_primary_cox <- function(dat, org, model_name, exposure_terms) {
  vars_needed <- c("followup_days", "adverse_event", "age", "sex", "race", "index_year_centered", "organ_score", "listing_center", exposure_terms)
  model_dat <- dat %>%
    filter(complete.cases(across(all_of(vars_needed))), followup_days > 0) %>%
    droplevels()
  if (nrow(model_dat) == 0L || sum(model_dat$adverse_event) < 100L) return(tibble())
  log_msg("Primary Cox ", org, " | ", model_name, " n=", nrow(model_dat), " adverse=", sum(model_dat$adverse_event))
  rhs <- c(exposure_terms, "age", "sex", "race", "index_year_centered", "organ_score", "listing_center", "cluster(PERS_ID)")
  form <- as.formula(paste("Surv(followup_days, adverse_event) ~", paste(rhs, collapse = " + ")))
  fit <- coxph(form, data = model_dat, ties = "efron", robust = TRUE, x = FALSE, y = FALSE)
  tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(term %in% exposure_terms) %>%
    transmute(
      organ = org,
      organ_label = recode(org, !!!organ_labels),
      model = model_name,
      endpoint = "death_or_deterioration_delist_cause_specific",
      term,
      exposure = recode(term, !!!term_labels),
      n = nrow(model_dat),
      adverse_events = sum(model_dat$adverse_event),
      competing_transplant_or_improvement = sum(model_dat$event_type == "transplant_or_improvement"),
      other_exits = sum(model_dat$event_type == "other_exit"),
      hazard_ratio = estimate,
      conf_low = conf.low,
      conf_high = conf.high,
      p_value = p.value,
      concordance = unname(summary(fit)$concordance[1]),
      variance_estimator = "robust_patient_clustered",
      center_adjustment = "listing_center_fixed_effect_dummies"
    )
}

log_msg("Fitting organ-specific formal primary Cox models")
model_result_dir <- file.path(out_dir, "primary_cox_model_results")
dir.create(model_result_dir, recursive = TRUE, showWarnings = FALSE)
cox_result_paths <- character()
for (org in target_organs) {
  dat <- analysis_dat %>% filter(WL_ORG == org)
  for (model_name in names(exposure_specs)) {
    result_path <- file.path(model_result_dir, paste0(tolower(org), "_", model_name, ".csv"))
    cox_result_paths <- c(cox_result_paths, result_path)
    if (file.exists(result_path)) {
      log_msg("Skipping existing primary Cox result ", org, " | ", model_name)
      next
    }
    result <- fit_primary_cox(dat, org, model_name, exposure_specs[[model_name]])
    write_csv(result, result_path)
    rm(result)
    invisible(gc())
  }
}

cox_results <- bind_rows(lapply(cox_result_paths[file.exists(cox_result_paths)], read_csv, show_col_types = FALSE))

write_csv(cox_results, file.path(out_dir, "formal_primary_cox_exposure_results.csv"))

write_csv(
  analysis_dat %>%
    group_by(WL_ORG, organ_score_name) %>%
    summarise(
      n = n(),
      people = n_distinct(PERS_ID),
      centers = n_distinct(listing_center),
      adverse_events = sum(adverse_event),
      transplant_or_improvement = sum(transplant_or_improvement),
      other_exit = sum(other_exit),
      median_followup_days = median(followup_days),
      complete_heat_humidity = sum(complete.cases(tmax_waitlist_5c, rmax_waitlist_10pct)),
      complete_pm25 = sum(!is.na(pm25_waitlist_5ug)),
      complete_o3 = sum(!is.na(o3_waitlist_10ppb)),
      complete_no2 = sum(!is.na(no2_waitlist_10unit)),
      mean_tmax_c = mean(tmax_c_waitlist_mean_exact_daily, na.rm = TRUE),
      mean_rmax_pct = mean(rmax_pct_waitlist_mean_exact_daily, na.rm = TRUE),
      mean_pm25 = mean(pm25_waitlist_ug_m3, na.rm = TRUE),
      mean_o3 = mean(o3_waitlist_ppb, na.rm = TRUE),
      mean_no2 = mean(no2_waitlist, na.rm = TRUE),
      .groups = "drop"
    ),
  file.path(out_dir, "formal_primary_cox_cohort_summary.csv")
)

log_msg("Wrote formal primary Cox outputs to ", out_dir)
