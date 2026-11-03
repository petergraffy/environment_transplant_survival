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
gridmet_dir <- file.path("output", "waitlist_environment_exposure", "gridmet_zcta_monthly_cache")
pollution_dir <- file.path("data", "release", "air_pollution_zcta_parquet")
in_dir <- file.path("output", "formal_waitlist_environment")
out_dir <- file.path("output", "formal_waitlist_environment_timevarying")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

analysis_start_year <- as.integer(Sys.getenv("ENV_ANALYSIS_START_YEAR", "2006"))
analysis_end_year <- as.integer(Sys.getenv("ENV_ANALYSIS_END_YEAR", "2023"))
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

flag_yes <- function(x) {
  y <- str_to_upper(str_trim(as.character(x)))
  as.integer(y %in% c("1", "Y", "YES", "TRUE", "T"))
}

safe_log1p <- function(x) {
  log(pmax(as.numeric(x), 0) + 1)
}

median_impute_num <- function(x) {
  x <- as.numeric(x)
  med <- median(x, na.rm = TRUE)
  if (!is.finite(med)) med <- NA_real_
  fifelse(is.na(x), med, x)
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

make_month_intervals <- function(interval_dt) {
  intervals <- interval_dt[, .(interval_id, waitlist_row_id, candidate_zip, interval_start_date, interval_end_date)]
  intervals[, month_seq := lapply(seq_len(.N), function(i) {
    start_month <- as.Date(format(interval_start_date[i], "%Y-%m-01"))
    end_month <- as.Date(format(interval_end_date[i], "%Y-%m-01"))
    if (is.na(start_month) || is.na(end_month) || end_month < start_month) return(as.Date(character()))
    seq(start_month, end_month, by = "month")
  })]
  intervals <- intervals[lengths(month_seq) > 0L]
  interval_lengths <- lengths(intervals$month_seq)
  out <- intervals[rep(seq_len(nrow(intervals)), interval_lengths)]
  out[, month_start := as.Date(unlist(intervals$month_seq, use.names = FALSE), origin = "1970-01-01")]
  out[, month_seq := NULL]
  out[, month_end := as.Date(format(month_start + 35L, "%Y-%m-01")) - 1L]
  out[, month_interval_start := pmax(interval_start_date, month_start)]
  out[, month_interval_end := pmin(interval_end_date, month_end)]
  out <- out[month_interval_end >= month_interval_start]
  out[, `:=`(
    year = as.integer(format(month_start, "%Y")),
    month = as.integer(format(month_start, "%m")),
    interval_days = as.integer(month_interval_end - month_interval_start) + 1L
  )]
  out
}

make_year_intervals <- function(interval_dt) {
  intervals <- interval_dt[, .(interval_id, waitlist_row_id, candidate_zip, interval_start_date, interval_end_date)]
  intervals[, year_seq := lapply(seq_len(.N), function(i) {
    start_year <- as.integer(format(interval_start_date[i], "%Y"))
    end_year <- as.integer(format(interval_end_date[i], "%Y"))
    if (is.na(start_year) || is.na(end_year) || end_year < start_year) return(integer())
    seq.int(max(start_year, analysis_start_year), min(end_year, analysis_end_year))
  })]
  intervals <- intervals[lengths(year_seq) > 0L]
  interval_lengths <- lengths(intervals$year_seq)
  out <- intervals[rep(seq_len(nrow(intervals)), interval_lengths)]
  out[, year := unlist(intervals$year_seq, use.names = FALSE)]
  out[, year_seq := NULL]
  out[, year_start := as.Date(paste0(year, "-01-01"))]
  out[, year_end := as.Date(paste0(year, "-12-31"))]
  out[, year_interval_start := pmax(interval_start_date, year_start)]
  out[, year_interval_end := pmin(interval_end_date, year_end)]
  out <- out[year_interval_end >= year_interval_start]
  out[, interval_days := as.integer(year_interval_end - year_interval_start) + 1L]
  out
}

candidate_files <- tribble(
  ~candidate_group, ~path,
  "kidney_pancreas", file.path(pubsaf_dir, "cand_kipa.sas7bdat"),
  "liver_intestine", file.path(pubsaf_dir, "cand_liin.sas7bdat"),
  "thoracic", file.path(pubsaf_dir, "cand_thor.sas7bdat")
)

stathist_files <- tribble(
  ~candidate_group, ~path,
  "kidney_pancreas", file.path(pubsaf_dir, "stathist_kipa.sas7bdat"),
  "liver_intestine", file.path(pubsaf_dir, "stathist_liin.sas7bdat"),
  "thoracic", file.path(pubsaf_dir, "stathist_thor.sas7bdat")
)

log_msg("Reading formal static analysis dataset")
base <- read_csv(file.path(in_dir, "formal_waitlist_environment_analysis_dataset.csv.gz"), show_col_types = FALSE) %>%
  filter(WL_ORG %in% target_organs) %>%
  mutate(
    event_type = factor(event_type, levels = c("censor", "adverse", "transplant_or_improvement", "other_exit")),
    sex = factor(sex),
    race = factor(race),
    listing_center = factor(listing_center),
    index_date = as.Date(index_date),
    observed_end_date = as.Date(observed_end_date)
  )

log_msg("Reconstructing candidate-level dialysis inputs")
candidate_extra <- candidate_files %>%
  rowwise() %>%
  mutate(data = list(read_sas(
    path,
    col_select = any_of(c("PERS_ID", "PX_ID", "WL_ORG", "CAN_ON_DIAL", "CAN_DIAL", "CAN_DIAL_DT"))
  ) %>% mutate(candidate_group = candidate_group))) %>%
  ungroup() %>%
  select(data) %>%
  unnest(data) %>%
  filter(WL_ORG %in% target_organs) %>%
  mutate(waitlist_row_id = row_number()) %>%
  transmute(
    waitlist_row_id,
    PERS_ID,
    PX_ID,
    WL_ORG,
    kidney_on_dialysis = flag_yes(coalesce(as.character(CAN_ON_DIAL), as.character(CAN_DIAL))),
    CAN_DIAL_DT
  )

base <- base %>%
  left_join(candidate_extra %>% select(waitlist_row_id, kidney_on_dialysis, CAN_DIAL_DT), by = "waitlist_row_id")

log_msg("Reading status-history intervals")
stathist <- stathist_files %>%
  rowwise() %>%
  mutate(data = list(read_sas(
    path,
    col_select = any_of(c(
      "PX_ID", "WL_ORG", "CANHX_BEGIN_DT", "CANHX_END_DT", "CANHX_STAT_CD",
      "CANHX_SRTR_LAB_MELD", "CANHX_OPTN_LAB_MELD"
    ))
  ))) %>%
  ungroup() %>%
  select(data) %>%
  unnest(data) %>%
  filter(WL_ORG %in% target_organs, !is.na(PX_ID), !is.na(CANHX_BEGIN_DT)) %>%
  as.data.table()

base_dt <- as.data.table(base)
setkey(base_dt, PX_ID, WL_ORG)
setkey(stathist, PX_ID, WL_ORG)

tv <- stathist[base_dt, allow.cartesian = TRUE, nomatch = 0L]
tv[, hx_end := fifelse(is.na(CANHX_END_DT), observed_end_date, CANHX_END_DT)]
tv[, interval_start_date := pmax(CANHX_BEGIN_DT, index_date)]
tv[, interval_end_date := pmin(hx_end, observed_end_date)]
tv <- tv[!is.na(interval_start_date) & !is.na(interval_end_date) & interval_end_date >= interval_start_date]
tv[, `:=`(
  tstart = as.numeric(interval_start_date - index_date),
  tstop = as.numeric(interval_end_date - index_date)
)]
tv[, tstop := pmax(tstop, tstart + 0.5)]
tv <- tv[tstop <= followup_days + 0.5]
tv[, tv_adverse_event := as.integer(adverse_event == 1L & interval_end_date >= observed_end_date)]
tv[, interval_id := .I]
tv[, current_inactive := as.integer(!is.na(CANHX_STAT_CD) & CANHX_STAT_CD %% 1000 == 999)]
tv[, current_status_unknown := as.integer(is.na(CANHX_STAT_CD) | CANHX_STAT_CD == 0)]

log_msg("Adding interval-updated organ score where available")
heart_updates <- read_sas(
  file.path(pubsaf_dir, "statjust_hr1a.sas7bdat"),
  col_select = any_of(c(
    "PX_ID", "WL_ORG", "CANHX_CHG_DT", "CANHX_LVAD_TYPE", "CANHX_RVAD_TYPE",
    "CANHX_ECMO", "CANHX_LAB_SERUM_CREAT", "CANHX_LAB_BILI",
    "CANHX_LAB_ALBUMIN", "CANHX_LAB_SODIUM", "CANHX_LAB_BNP"
  ))
) %>%
  filter(WL_ORG == "HR", !is.na(PX_ID), !is.na(CANHX_CHG_DT)) %>%
  as.data.table()

if (nrow(heart_updates) > 0L) {
  heart_updates[, `:=`(
    hr_albumin_imputed = median_impute_num(CANHX_LAB_ALBUMIN),
    hr_bilirubin_imputed = median_impute_num(CANHX_LAB_BILI),
    hr_creatinine_imputed = median_impute_num(CANHX_LAB_SERUM_CREAT),
    hr_sodium_imputed = median_impute_num(CANHX_LAB_SODIUM),
    hr_bnp_imputed = median_impute_num(CANHX_LAB_BNP),
    hr_short_mcs = as.integer(flag_yes(CANHX_ECMO) == 1L | flag_yes(CANHX_RVAD_TYPE) == 1L),
    hr_durable_lvad = flag_yes(CANHX_LVAD_TYPE)
  )]
  heart_updates[, hr_egfr := calc_egfr_2021(hr_creatinine_imputed, 55, "M")]
  heart_updates[, heart_score_update := -0.656 * hr_albumin_imputed +
    0.617 * safe_log1p(hr_bilirubin_imputed) -
    0.012 * hr_egfr -
    0.077 * hr_sodium_imputed -
    0.377 * hr_durable_lvad +
    1.092 * hr_short_mcs +
    0.433 * safe_log1p(hr_bnp_imputed)]
  setkey(heart_updates, PX_ID, CANHX_CHG_DT)
  hr_intervals <- tv[WL_ORG == "HR", .(interval_id, PX_ID, interval_start_date)]
  setkey(hr_intervals, PX_ID, interval_start_date)
  hr_join <- heart_updates[hr_intervals, on = c("PX_ID", "CANHX_CHG_DT" = "interval_start_date"), roll = Inf]
  tv[hr_join, heart_score_update := i.heart_score_update, on = "interval_id"]
}

tv[, age_interval := age + tstart / 365.25]
tv[, dialysis_years_interval := as.numeric(interval_start_date - CAN_DIAL_DT) / 365.25]
tv[is.na(dialysis_years_interval) | dialysis_years_interval < 0, dialysis_years_interval := 0]
tv[, kidney_score_update := age_interval / 10 + 0.8 * kidney_on_dialysis + 0.3 * log1p(dialysis_years_interval)]
tv[, liver_score_update := fifelse(!is.na(CANHX_SRTR_LAB_MELD), CANHX_SRTR_LAB_MELD, CANHX_OPTN_LAB_MELD)]
tv[, organ_score_tv := fcase(
  WL_ORG == "HR" & !is.na(heart_score_update), heart_score_update,
  WL_ORG == "KI" & !is.na(kidney_score_update), kidney_score_update,
  WL_ORG == "LI" & !is.na(liver_score_update), liver_score_update,
  default = organ_score
)]

log_msg("Aggregating interval exposures with month/day weights")
month_intervals <- make_month_intervals(tv)
year_intervals <- make_year_intervals(tv)

gridmet_month_parts <- list()
for (yr in analysis_start_year:analysis_end_year) {
  cache_path <- file.path(gridmet_dir, sprintf("gridmet_zcta_monthly_%04d.csv.gz", yr))
  if (!file.exists(cache_path)) next
  intervals_y <- month_intervals[year == yr]
  if (nrow(intervals_y) == 0L) next
  gm <- read_csv(cache_path, show_col_types = FALSE, col_types = cols(zip = col_character(), .default = col_double())) %>%
    transmute(
      zip = clean_zip(zip),
      year,
      month,
      tmax_c,
      rmax_pct
    ) %>%
    as.data.table()
  joined <- gm[intervals_y, on = c("zip" = "candidate_zip", "year", "month")]
  gridmet_month_parts[[as.character(yr)]] <- joined[!is.na(tmax_c) & !is.na(rmax_pct), .(
    tmax_days = sum(interval_days),
    rmax_days = sum(interval_days),
    tmax_interval_c = sum(tmax_c * interval_days) / sum(interval_days),
    rmax_interval_pct = sum(rmax_pct * interval_days) / sum(interval_days)
  ), by = interval_id]
}
gridmet_interval <- rbindlist(gridmet_month_parts, use.names = TRUE, fill = TRUE)[
  ,
  .(
    tmax_days = sum(tmax_days),
    rmax_days = sum(rmax_days),
    tmax_interval_c = weighted.mean(tmax_interval_c, tmax_days),
    rmax_interval_pct = weighted.mean(rmax_interval_pct, rmax_days)
  ),
  by = interval_id
]

pm25 <- read_parquet(file.path(pollution_dir, "air_pollution_zcta_pm25_monthly_2005_2023.parquet")) %>%
  transmute(zip = clean_zip(zip), year, month, pm25_ug_m3) %>%
  as.data.table()
o3 <- read_parquet(file.path(pollution_dir, "air_pollution_zcta_o3_monthly_2005_2023.parquet")) %>%
  transmute(zip = clean_zip(zip), year, month, o3_ppb) %>%
  as.data.table()
no2 <- read_parquet(file.path(pollution_dir, "air_pollution_zcta_no2_annual_2005_2025.parquet")) %>%
  transmute(zip = clean_zip(zip), year, no2) %>%
  as.data.table()

pm25_interval <- pm25[month_intervals, on = c("zip" = "candidate_zip", "year", "month")][!is.na(pm25_ug_m3), .(
  pm25_interval_ug_m3 = sum(pm25_ug_m3 * interval_days) / sum(interval_days)
), by = interval_id]
o3_interval <- o3[month_intervals, on = c("zip" = "candidate_zip", "year", "month")][!is.na(o3_ppb), .(
  o3_interval_ppb = sum(o3_ppb * interval_days) / sum(interval_days)
), by = interval_id]
no2_interval <- no2[year_intervals, on = c("zip" = "candidate_zip", "year")][!is.na(no2), .(
  no2_interval = sum(no2 * interval_days) / sum(interval_days)
), by = interval_id]

setkey(tv, interval_id)
for (obj in list(gridmet_interval, pm25_interval, o3_interval, no2_interval)) {
  setkey(obj, interval_id)
  tv <- obj[tv]
}

tv[, `:=`(
  tmax_interval_5c = tmax_interval_c / 5,
  rmax_interval_10pct = rmax_interval_pct / 10,
  pm25_interval_5ug = pm25_interval_ug_m3 / 5,
  o3_interval_10ppb = o3_interval_ppb / 10,
  no2_interval_10unit = no2_interval / 10
)]

write_csv(
  as_tibble(tv[, .(
    interval_id, waitlist_row_id, PERS_ID, PX_ID, WL_ORG, candidate_zip,
    interval_start_date, interval_end_date, tstart, tstop, tv_adverse_event,
    age_interval, sex, race, index_year_centered, listing_center, organ_score_tv,
    current_inactive, current_status_unknown,
    tmax_interval_c, rmax_interval_pct, pm25_interval_ug_m3, o3_interval_ppb, no2_interval
  )]),
  file.path(out_dir, "timevarying_interval_analysis_dataset.csv.gz")
)

if (Sys.getenv("BUILD_INTERVALS_ONLY", "1") == "1") {
  log_msg("Wrote time-varying interval dataset to ", out_dir)
  log_msg("Set BUILD_INTERVALS_ONLY=0 to run the legacy all-in-one modeling branch; preferred fitters are code/33 and code/34.")
  quit(save = "no", status = 0)
}

exposure_specs <- list(
  heat_humidity = c("tmax_interval_5c", "rmax_interval_10pct"),
  pm25 = "pm25_interval_5ug",
  o3 = "o3_interval_10ppb",
  no2 = "no2_interval_10unit"
)

term_labels <- c(
  tmax_interval_5c = "Interval Tmax per 5 C",
  rmax_interval_10pct = "Interval maximum relative humidity per 10 pct",
  pm25_interval_5ug = "Interval PM2.5 per 5 ug/m3",
  o3_interval_10ppb = "Interval ozone per 10 ppb",
  no2_interval_10unit = "Interval NO2 per 10 units"
)

fit_tv_cox <- function(dat, org, model_name, exposure_terms) {
  vars_needed <- c("tstart", "tstop", "tv_adverse_event", "age_interval", "sex", "race", "index_year_centered", "organ_score_tv", "listing_center", exposure_terms)
  model_dat <- dat[complete.cases(dat[, ..vars_needed]) & tstop > tstart]
  if (nrow(model_dat) == 0L || sum(model_dat$tv_adverse_event) < 100L) return(tibble())
  log_msg("Time-varying Cox ", org, " | ", model_name, " intervals=", nrow(model_dat), " adverse=", sum(model_dat$tv_adverse_event))
  rhs <- c(exposure_terms, "age_interval", "sex", "race", "index_year_centered", "organ_score_tv", "strata(listing_center)", "cluster(PERS_ID)")
  form <- as.formula(paste("Surv(tstart, tstop, tv_adverse_event) ~", paste(rhs, collapse = " + ")))
  fit <- coxph(form, data = model_dat, ties = "efron", robust = TRUE, x = FALSE, y = FALSE)
  tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(term %in% exposure_terms) %>%
    transmute(
      organ = org,
      organ_label = recode(org, !!!organ_labels),
      model = model_name,
      endpoint = "time_varying_death_or_deterioration_delist_cause_specific",
      term,
      exposure = recode(term, !!!term_labels),
      intervals = nrow(model_dat),
      waitlist_rows = uniqueN(model_dat$waitlist_row_id),
      people = uniqueN(model_dat$PERS_ID),
      adverse_events = sum(model_dat$tv_adverse_event),
      hazard_ratio = estimate,
      conf_low = conf.low,
      conf_high = conf.high,
      p_value = p.value,
      concordance = unname(summary(fit)$concordance[1]),
      variance_estimator = "robust_patient_clustered",
      center_adjustment = "stratified_listing_center_baseline_hazard"
    )
}

model_result_dir <- file.path(out_dir, "timevarying_cox_model_results")
dir.create(model_result_dir, recursive = TRUE, showWarnings = FALSE)
tv_result_paths <- character()
for (org in target_organs) {
  dat <- tv[WL_ORG == org]
  for (model_name in names(exposure_specs)) {
    result_path <- file.path(model_result_dir, paste0(tolower(org), "_", model_name, ".csv"))
    tv_result_paths <- c(tv_result_paths, result_path)
    if (file.exists(result_path)) {
      log_msg("Skipping existing time-varying Cox result ", org, " | ", model_name)
      next
    }
    result <- fit_tv_cox(dat, org, model_name, exposure_specs[[model_name]])
    write_csv(result, result_path)
    rm(result)
    invisible(gc())
  }
}

tv_results <- bind_rows(lapply(tv_result_paths[file.exists(tv_result_paths)], read_csv, show_col_types = FALSE))
write_csv(tv_results, file.path(out_dir, "timevarying_cox_exposure_results.csv"))

log_msg("Fitting reduced Fine-Gray sensitivity with month/day-weighted Tmax and Rmax")
fg_base <- tv[
  ,
  .(
    followup_days = max(tstop, na.rm = TRUE),
    event_type = first(event_type),
    age = first(age),
    sex = first(sex),
    race = first(race),
    index_year_centered = first(index_year_centered),
    organ_score = first(organ_score),
    listing_center = first(listing_center),
    tmax_month_weighted_5c = weighted.mean(tmax_interval_5c, pmax(tstop - tstart, 0.5), na.rm = TRUE),
    rmax_month_weighted_10pct = weighted.mean(rmax_interval_10pct, pmax(tstop - tstart, 0.5), na.rm = TRUE)
  ),
  by = .(waitlist_row_id, PERS_ID, WL_ORG)
]

fit_reduced_fg <- function(dat, org) {
  model_dat <- dat[
    WL_ORG == org &
      complete.cases(followup_days, event_type, age, sex, race, index_year_centered, organ_score, tmax_month_weighted_5c, rmax_month_weighted_10pct)
  ]
  if (nrow(model_dat) == 0L || sum(model_dat$event_type == "adverse") < 100L) return(tibble())
  log_msg("Reduced Fine-Gray heat/humidity ", org, " n=", nrow(model_dat), " adverse=", sum(model_dat$event_type == "adverse"))
  fstatus <- fcase(
    model_dat$event_type == "adverse", 1L,
    model_dat$event_type %in% c("transplant_or_improvement", "other_exit"), 2L,
    default = 0L
  )
  mm <- model.matrix(
    ~ tmax_month_weighted_5c + rmax_month_weighted_10pct + age + sex + race + index_year_centered + organ_score,
    data = model_dat
  )
  mm <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]
  fit <- crr(model_dat$followup_days, fstatus, cov1 = mm, failcode = 1, cencode = 0)
  coef <- fit$coef
  var <- diag(fit$var)
  keep <- intersect(c("tmax_month_weighted_5c", "rmax_month_weighted_10pct"), names(coef))
  tibble(
    organ = org,
    organ_label = recode(org, !!!organ_labels),
    model = "reduced_finegray_heat_humidity_month_weighted",
    endpoint = "death_or_deterioration_delist_subdistribution",
    term = keep,
    exposure = recode(
      keep,
      tmax_month_weighted_5c = "Month/day-weighted Tmax per 5 C",
      rmax_month_weighted_10pct = "Month/day-weighted maximum relative humidity per 10 pct"
    ),
    n = nrow(model_dat),
    adverse_events = sum(model_dat$event_type == "adverse"),
    subdistribution_hazard_ratio = exp(coef[keep]),
    conf_low = exp(coef[keep] - 1.96 * sqrt(var[match(keep, names(coef))])),
    conf_high = exp(coef[keep] + 1.96 * sqrt(var[match(keep, names(coef))])),
    p_value = 2 * pnorm(abs(coef[keep] / sqrt(var[match(keep, names(coef))])), lower.tail = FALSE),
    center_adjustment = "not_included_reduced_sensitivity",
    variance_estimator = "cmprsk_crr_model_based"
  )
}

fg_results <- bind_rows(lapply(target_organs, function(org) fit_reduced_fg(fg_base, org)))
write_csv(fg_results, file.path(out_dir, "reduced_finegray_heat_humidity_results.csv"))

write_csv(
  as_tibble(tv[, .(
    intervals = .N,
    waitlist_rows = uniqueN(waitlist_row_id),
    people = uniqueN(PERS_ID),
    adverse_events = sum(tv_adverse_event),
    median_intervals_per_waitlist = median(.N),
    complete_heat_humidity_intervals = sum(complete.cases(tmax_interval_5c, rmax_interval_10pct)),
    complete_pm25_intervals = sum(!is.na(pm25_interval_5ug)),
    complete_o3_intervals = sum(!is.na(o3_interval_10ppb)),
    complete_no2_intervals = sum(!is.na(no2_interval_10unit))
  ), by = WL_ORG]),
  file.path(out_dir, "timevarying_interval_summary.csv")
)

log_msg("Wrote time-varying Cox and reduced Fine-Gray outputs to ", out_dir)
