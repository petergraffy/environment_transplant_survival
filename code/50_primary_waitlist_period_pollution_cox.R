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
saf_paths <- get_saf_paths(release = "q1_2026")
assert_saf_files(saf_paths, include_stathist = TRUE)

pubsaf_dir <- saf_paths$pubsaf_dir
pollution_dir <- file.path("data", "release", "air_pollution_zcta_parquet")
out_dir <- file.path("output", "primary_waitlist_period_pollution_cox")
model_result_dir <- file.path(out_dir, "model_results")
dir.create(model_result_dir, recursive = TRUE, showWarnings = FALSE)

analysis_start_year <- as.integer(Sys.getenv("ENV_ANALYSIS_START_YEAR", "2005"))
pollution_end_year <- c(pm25 = 2023L, o3 = 2023L, no2 = 2025L)
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

durable_vad_brand_codes <- c(
  202L, 205L, 206L, 207L, 208L, 209L, 210L, 212L, 213L, 214L,
  223L, 224L, 233L, 236L, 239L, 240L, 312L, 313L, 314L, 315L,
  316L, 319L, 322L, 327L, 330L, 333L, 334L
)

flag_durable_vad_brand <- function(...) {
  vals <- list(...)
  Reduce(`|`, lapply(vals, function(x) suppressWarnings(as.integer(as.character(x))) %in% durable_vad_brand_codes)) %>%
    as.integer()
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

make_year_intervals <- function(cohort_dt, exposure_end_year) {
  intervals <- cohort_dt[, .(waitlist_row_id, candidate_zip, index_date, observed_end_date)]
  intervals[, interval_years := lapply(seq_len(.N), function(i) {
    start_year <- max(as.integer(format(index_date[i], "%Y")), analysis_start_year)
    stop_year <- min(as.integer(format(observed_end_date[i], "%Y")), exposure_end_year)
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

candidate_files <- tribble(
  ~candidate_group, ~path,
  "kidney_pancreas", file.path(pubsaf_dir, "cand_kipa.sas7bdat"),
  "liver_intestine", file.path(pubsaf_dir, "cand_liin.sas7bdat"),
  "thoracic", file.path(pubsaf_dir, "cand_thor.sas7bdat")
)

candidate_cols <- c(
  "PERS_ID", "PX_ID", "WL_ORG", "CAN_LISTING_DT", "CAN_ACTIVATE_DT", "CAN_SOURCE",
  "CAN_REM_DT", "CAN_REM_CD", "CAN_DEATH_DT", "CAN_ENDWLFU", "CAN_AGE_AT_LISTING",
  "CAN_GENDER", "CAN_RACE_SRTR", "CAN_LISTING_CTR_CD", "CAN_ON_DIAL",
  "CAN_DIAL", "CAN_DIAL_DT", "CAN_DIAB", "CAN_DIAB_TY", "CAN_MOST_RECENT_CREAT", "CAN_INIT_SRTR_LAB_MELD",
  "CAN_LAST_SRTR_LAB_MELD", "CAN_TOT_BILI", "CAN_TOT_ALBUMIN",
  "CAN_LAST_SERUM_SODIUM", "CAN_FVC", "CAN_FEV1", "CAN_PCO2", "CAN_AT_REST_O2",
  "CAN_SIX_MIN_WALK", "CAN_PULM_ART_MEAN", "CAN_CARDIAC_OUTPUT",
  "CAN_VENTILATOR", "CAN_ON_VENTILATOR", "CAN_ECMO", "CAN_VAD1",
  "CAN_VAD2", "CAN_CORTICOST_DEPND"
)

log_msg("Reading Q1 2026 SAF candidate and ZIP inputs from ", saf_paths$saf_dir)
candidate_zip <- read_sas(saf_paths$canzip_file) %>%
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
    "PX_ID", "WL_ORG", "CANHX_CHG_DT", "CANHX_RVAD_TYPE",
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
    hr_creatinine = CANHX_LAB_SERUM_CREAT,
    hr_bilirubin = CANHX_LAB_BILI,
    hr_albumin = CANHX_LAB_ALBUMIN,
    hr_sodium = CANHX_LAB_SODIUM,
    hr_bnp = CANHX_LAB_BNP
  )

log_msg("Constructing waitlist cohort")
max_end_date <- as.Date(sprintf("%04d-12-31", max(pollution_end_year)))
cohort <- candidate_data %>%
  filter(WL_ORG %in% target_organs) %>%
  left_join(candidate_zip, by = c("PERS_ID", "PX_ID")) %>%
  left_join(heart_baseline, by = "PX_ID") %>%
  mutate(
    waitlist_row_id = row_number(),
    index_date = coalesce(CAN_ACTIVATE_DT, CAN_LISTING_DT),
    index_year = as.integer(format(index_date, "%Y")),
    raw_event_date = coalesce(CAN_REM_DT, CAN_DEATH_DT, CAN_ENDWLFU),
    observed_end_date = if_else(is.na(raw_event_date) | raw_event_date > max_end_date, max_end_date, raw_event_date),
    event_observed = !is.na(raw_event_date) & raw_event_date <= max_end_date,
    adverse_event = as.integer(event_observed & CAN_REM_CD %in% c(8, 13)),
    transplant_or_improvement = as.integer(event_observed & CAN_REM_CD %in% c(4, 12, 14, 15, 18, 19, 21, 22, 23)),
    other_exit = as.integer(event_observed & adverse_event == 0L & transplant_or_improvement == 0L & !is.na(CAN_REM_CD)),
    event_type = factor(
      case_when(
        adverse_event == 1L ~ "death_or_deterioration_delist",
        transplant_or_improvement == 1L ~ "transplant_or_improvement",
        other_exit == 1L ~ "other_exit",
        TRUE ~ "administrative_censor"
      ),
      levels = c("administrative_censor", "death_or_deterioration_delist", "transplant_or_improvement", "other_exit")
    ),
    age = as.numeric(CAN_AGE_AT_LISTING),
    sex = as_factor_missing(if_else(CAN_GENDER %in% c("M", "F"), CAN_GENDER, NA_character_)),
    race = as_factor_missing(CAN_RACE_SRTR),
    listing_center = as_factor_missing(CAN_LISTING_CTR_CD),
    listing_year = factor(index_year),
    kidney_on_dialysis = flag_yes(coalesce(as.character(CAN_ON_DIAL), as.character(CAN_DIAL))),
    kidney_diabetes = as.integer(
      CAN_DIAB %in% c(2, 3, 4) |
        CAN_DIAB_TY %in% c(2, 3, 4, 5)
    ),
    dialysis_years_raw = as.numeric(index_date - CAN_DIAL_DT) / 365.25,
    dialysis_years_raw = if_else(is.na(dialysis_years_raw) | dialysis_years_raw < 0, NA_real_, dialysis_years_raw),
    kidney_no_dialysis_time = as.integer(kidney_on_dialysis != 1L | is.na(dialysis_years_raw) | dialysis_years_raw == 0),
    kidney_dialysis_years = if_else(kidney_no_dialysis_time == 1L, 0, dialysis_years_raw),
    dialysis_years = kidney_dialysis_years,
    liver_meld = coalesce(CAN_INIT_SRTR_LAB_MELD, CAN_LAST_SRTR_LAB_MELD),
    hr_creatinine_for_score = coalesce(hr_creatinine, CAN_MOST_RECENT_CREAT),
    hr_albumin = coalesce(hr_albumin, CAN_TOT_ALBUMIN),
    hr_bilirubin = coalesce(hr_bilirubin, CAN_TOT_BILI),
    hr_sodium = coalesce(hr_sodium, CAN_LAST_SERUM_SODIUM),
    hr_short_mcs = coalesce(hr_short_mcs, flag_yes(CAN_ECMO), 0L),
    hr_durable_lvad = if_else(WL_ORG == "HR", flag_durable_vad_brand(CAN_VAD1, CAN_VAD2), 0L),
    lung_fev1_score_input = as.numeric(coalesce(CAN_FEV1, 0)),
    lung_fvc_score_input = as.numeric(coalesce(CAN_FVC, 0)),
    lung_pco2_score_input = as.numeric(coalesce(CAN_PCO2, 0)),
    lung_resting_o2_score_input = as.numeric(coalesce(CAN_AT_REST_O2, 0)),
    lung_six_min_walk_score_input = as.numeric(coalesce(CAN_SIX_MIN_WALK, 0)),
    lung_pulm_art_mean_score_input = as.numeric(coalesce(CAN_PULM_ART_MEAN, 0)),
    lung_cardiac_output_score_input = as.numeric(coalesce(CAN_CARDIAC_OUTPUT, 0)),
    lung_ventilator_score_input = flag_yes(coalesce(as.character(CAN_VENTILATOR), as.character(CAN_ON_VENTILATOR))),
    lung_ecmo_score_input = flag_yes(CAN_ECMO),
    lung_corticosteroid_score_input = flag_yes(CAN_CORTICOST_DEPND)
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
      lung_fev1_score_input * -0.01,
      lung_fvc_score_input * -0.005,
      lung_pco2_score_input * 0.02,
      lung_resting_o2_score_input * 0.01,
      lung_six_min_walk_score_input * -0.002,
      lung_pulm_art_mean_score_input * 0.01,
      lung_cardiac_output_score_input * -0.05,
      lung_ventilator_score_input * 1.0,
      lung_ecmo_score_input * 1.5,
      lung_corticosteroid_score_input * 0.3
    ), na.rm = TRUE),
    organ_score = case_when(
      WL_ORG == "HR" ~ us_crs_proxy,
      WL_ORG == "KI" ~ NA_real_,
      WL_ORG == "LI" ~ liver_meld,
      WL_ORG == "LU" ~ lung_las_cas_component_score
    ),
    organ_score_name = case_when(
      WL_ORG == "HR" ~ "US-CRS proxy baseline",
      WL_ORG == "KI" ~ "Age + dialysis + diabetes baseline covariates",
      WL_ORG == "LI" ~ "Initial MELD, last if initial missing",
      WL_ORG == "LU" ~ "LAS/CAS component proxy baseline"
    )
  ) %>%
  filter(
    index_year >= analysis_start_year,
    index_year <= max(pollution_end_year),
    CAN_SOURCE %in% c("A", "R"),
    !is.na(index_date),
    !is.na(observed_end_date),
    observed_end_date >= index_date,
    !is.na(candidate_zip)
  )

cohort_dt <- as.data.table(cohort)

log_msg("Reading pollution release tables")
pm25 <- read_parquet(file.path(pollution_dir, "air_pollution_zcta_pm25_monthly_2005_2023.parquet")) %>%
  transmute(zip = clean_zip(zip), year = as.integer(year), pm25_ug_m3) %>%
  group_by(zip, year) %>%
  summarise(pm25_n_months = sum(!is.na(pm25_ug_m3)), pm25_ug_m3 = mean(pm25_ug_m3, na.rm = TRUE), .groups = "drop") %>%
  filter(pm25_n_months == 12L, is.finite(pm25_ug_m3)) %>%
  as.data.table()
o3 <- read_parquet(file.path(pollution_dir, "air_pollution_zcta_o3_monthly_2005_2023.parquet")) %>%
  transmute(zip = clean_zip(zip), year = as.integer(year), o3_ppb) %>%
  group_by(zip, year) %>%
  summarise(o3_n_months = sum(!is.na(o3_ppb)), o3_ppb = mean(o3_ppb, na.rm = TRUE), .groups = "drop") %>%
  filter(o3_n_months == 12L, is.finite(o3_ppb)) %>%
  as.data.table()
no2 <- read_parquet(file.path(pollution_dir, "air_pollution_zcta_no2_annual_2005_2025.parquet")) %>%
  transmute(zip = clean_zip(zip), year, no2) %>%
  as.data.table()

log_msg("Computing year/day-weighted waitlist-period PM2.5 and ozone means through 2023")
year_intervals_2023 <- make_year_intervals(cohort_dt, exposure_end_year = 2023L)
pm25_year <- pm25[year_intervals_2023, on = c("zip" = "candidate_zip", "year" = "interval_year")]
o3_year <- o3[year_intervals_2023, on = c("zip" = "candidate_zip", "year" = "interval_year")]
pm25_exp <- pm25_year[!is.na(pm25_ug_m3), .(
  pm25_days = sum(interval_days),
  pm25_waitlist_ug_m3 = sum(pm25_ug_m3 * interval_days) / sum(interval_days)
), by = waitlist_row_id]
o3_exp <- o3_year[!is.na(o3_ppb), .(
  o3_days = sum(interval_days),
  o3_waitlist_ppb = sum(o3_ppb * interval_days) / sum(interval_days)
), by = waitlist_row_id]

log_msg("Computing waitlist-period NO2 means through 2025")
year_intervals_2025 <- make_year_intervals(cohort_dt, exposure_end_year = 2025L)
no2_year <- no2[year_intervals_2025, on = c("zip" = "candidate_zip", "year" = "interval_year")]
no2_exp <- no2_year[!is.na(no2), .(
  no2_days = sum(interval_days),
  no2_waitlist = sum(no2 * interval_days) / sum(interval_days)
), by = waitlist_row_id]

analysis_dat <- cohort %>%
  left_join(as_tibble(pm25_exp), by = "waitlist_row_id") %>%
  left_join(as_tibble(o3_exp), by = "waitlist_row_id") %>%
  left_join(as_tibble(no2_exp), by = "waitlist_row_id") %>%
  mutate(
    followup_days_pm25 = as.numeric(pmin(observed_end_date, as.Date("2023-12-31")) - index_date),
    followup_days_o3 = followup_days_pm25,
    followup_days_no2 = as.numeric(pmin(observed_end_date, as.Date("2025-12-31")) - index_date),
    event_pm25 = as.integer(adverse_event == 1L & observed_end_date <= as.Date("2023-12-31")),
    event_o3 = event_pm25,
    event_no2 = adverse_event,
    pm25_waitlist_5ug = pm25_waitlist_ug_m3 / 5,
    o3_waitlist_10ppb = o3_waitlist_ppb / 10,
    no2_waitlist_10unit = no2_waitlist / 10
  )

write_csv(
  analysis_dat %>%
    select(
      waitlist_row_id, PERS_ID, PX_ID, WL_ORG, candidate_zip, index_date, observed_end_date,
      event_type, adverse_event, transplant_or_improvement, other_exit,
      age, sex, race, listing_year, organ_score, organ_score_name, listing_center,
      hr_albumin_imputed, hr_bilirubin_imputed, hr_egfr, hr_sodium_imputed,
      hr_durable_lvad, hr_short_mcs, hr_bnp_imputed,
      kidney_on_dialysis, kidney_no_dialysis_time, kidney_dialysis_years, kidney_diabetes, dialysis_years, liver_meld,
      lung_fev1_score_input, lung_fvc_score_input, lung_pco2_score_input,
      lung_resting_o2_score_input, lung_six_min_walk_score_input,
      lung_pulm_art_mean_score_input, lung_cardiac_output_score_input,
      lung_ventilator_score_input, lung_ecmo_score_input, lung_corticosteroid_score_input,
      followup_days_pm25, event_pm25, pm25_waitlist_ug_m3, pm25_days,
      followup_days_o3, event_o3, o3_waitlist_ppb, o3_days,
      followup_days_no2, event_no2, no2_waitlist, no2_days
    ),
  file.path(out_dir, "primary_waitlist_period_pollution_analysis_dataset.csv.gz")
)

exposure_specs <- list(
  pm25 = list(term = "pm25_waitlist_5ug", label = "PM2.5 waitlist-period mean per 5 ug/m3", followup = "followup_days_pm25", event = "event_pm25", end_year = 2023L),
  o3 = list(term = "o3_waitlist_10ppb", label = "Ozone waitlist-period mean per 10 ppb", followup = "followup_days_o3", event = "event_o3", end_year = 2023L),
  no2 = list(term = "no2_waitlist_10unit", label = "NO2 waitlist-period mean per 10-unit increment", followup = "followup_days_no2", event = "event_no2", end_year = 2025L)
)

fit_waitlist_cox <- function(dat, org, exposure_name, spec) {
  adjustment_terms <- c("age", "sex", "race", "listing_year", "strata(listing_center)")
  adjustment_vars <- c("age", "sex", "race", "listing_year", "listing_center")
  if (org == "KI") {
    adjustment_terms <- c("age", "sex", "race", "listing_year", "kidney_no_dialysis_time", "kidney_dialysis_years", "kidney_diabetes", "strata(listing_center)")
    adjustment_vars <- c("age", "sex", "race", "listing_year", "kidney_no_dialysis_time", "kidney_dialysis_years", "kidney_diabetes", "listing_center")
  } else {
    adjustment_terms <- c("age", "sex", "race", "listing_year", "organ_score", "strata(listing_center)")
    adjustment_vars <- c("age", "sex", "race", "listing_year", "organ_score", "listing_center")
  }
  vars_needed <- c(spec$followup, spec$event, adjustment_vars, spec$term)
  model_dat <- dat %>%
    filter(
      index_year <= spec$end_year,
      complete.cases(across(all_of(vars_needed))),
      .data[[spec$followup]] > 0
    ) %>%
    mutate(.followup_days = pmax(.data[[spec$followup]], 0.5), .event = .data[[spec$event]]) %>%
    droplevels()

  if (nrow(model_dat) == 0L || sum(model_dat$.event) < 50L || n_distinct(model_dat$listing_center) < 2L) {
    return(tibble())
  }

  log_msg("Waitlist-period pollution Cox ", org, " | ", exposure_name, " n=", nrow(model_dat), " adverse=", sum(model_dat$.event))
  form <- as.formula(paste(
    "Surv(.followup_days, .event) ~",
    paste(c(spec$term, adjustment_terms), collapse = " + ")
  ))
  fit <- coxph(form, data = model_dat, ties = "efron", x = FALSE, y = FALSE)

  tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(term == spec$term) %>%
    transmute(
      organ = org,
      organ_label = recode(org, !!!organ_labels),
      exposure = exposure_name,
      exposure_label = spec$label,
      exposure_window = "observed_waitlist_period_mean",
      exposure_data_end_year = spec$end_year,
      endpoint = "death_or_deterioration_delist_cause_specific",
      n = nrow(model_dat),
      people = n_distinct(model_dat$PERS_ID),
      centers = n_distinct(model_dat$listing_center),
      adverse_events = sum(model_dat$.event),
      hazard_ratio = estimate,
      conf_low = conf.low,
      conf_high = conf.high,
      p_value = p.value,
      concordance = unname(summary(fit)$concordance[1]),
      adjustment_set = if_else(
        org == "KI",
        "age + sex + race + listing_year + no_dialysis_time + dialysis_duration_years + diabetes + strata(listing_center)",
        "age + sex + race + listing_year + baseline_organ_score + strata(listing_center)"
      ),
      variance_estimator = "model_based",
      center_adjustment = "stratified_baseline_hazard_by_listing_center"
    )
}

log_msg("Fitting organ-specific primary waitlist-period pollution Cox models")
result_paths <- character()
for (org in target_organs) {
  org_dat <- analysis_dat %>% filter(WL_ORG == org)
  for (exposure_name in names(exposure_specs)) {
    result_path <- file.path(model_result_dir, paste0(tolower(org), "_", exposure_name, ".csv"))
    result_paths <- c(result_paths, result_path)
    result <- fit_waitlist_cox(org_dat, org, exposure_name, exposure_specs[[exposure_name]])
    write_csv(result, result_path)
    rm(result)
    invisible(gc())
  }
}

waitlist_results <- bind_rows(lapply(result_paths[file.exists(result_paths)], read_csv, show_col_types = FALSE))
write_csv(waitlist_results, file.path(out_dir, "primary_waitlist_period_pollution_cox_results.csv"))

write_csv(
  bind_rows(lapply(names(exposure_specs), function(exposure_name) {
    spec <- exposure_specs[[exposure_name]]
    analysis_dat %>%
      filter(index_year <= spec$end_year, !is.na(.data[[spec$term]]), .data[[spec$followup]] > 0) %>%
      group_by(WL_ORG, organ_score_name) %>%
      summarise(
        exposure = exposure_name,
        exposure_data_end_year = spec$end_year,
        n = n(),
        people = n_distinct(PERS_ID),
        centers = n_distinct(listing_center),
        adverse_events = sum(.data[[spec$event]], na.rm = TRUE),
        median_followup_days = median(.data[[spec$followup]], na.rm = TRUE),
        mean_exposure = mean(.data[[spec$term]], na.rm = TRUE),
        .groups = "drop"
      )
  })),
  file.path(out_dir, "primary_waitlist_period_pollution_cohort_summary.csv")
)

log_msg("Wrote primary waitlist-period pollution Cox outputs to ", out_dir)
