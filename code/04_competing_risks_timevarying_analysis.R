#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(broom)
  library(cmprsk)
  library(data.table)
  library(dplyr)
  library(haven)
  library(readr)
  library(splines)
  library(stringr)
  library(survival)
  library(tibble)
})

saf_paths_source <- file.path("code", "saf_paths.R")
if (!file.exists(saf_paths_source)) saf_paths_source <- "saf_paths.R"
source(saf_paths_source)
saf_paths <- get_saf_paths()
assert_saf_files(saf_paths, include_stathist = TRUE)
saf_dir <- saf_paths$saf_dir
pubsaf_dir <- saf_paths$pubsaf_dir
supp_dir <- saf_paths$supp_dir
pollution_dir <- "output/zip_pollution"
out_dir <- "output/competing_timevarying"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
message("Using SAF directory: ", saf_dir)

analysis_start_year <- 2005L
analysis_end_year <- 2023L

log_msg <- function(...) {
  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
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

clean_zip <- function(x) {
  z <- str_extract(as.character(x), "[0-9]+")
  z <- if_else(nchar(z) > 5L, substr(z, 1, 5), z)
  if_else(is.na(z), NA_character_, str_pad(z, 5, side = "left", pad = "0"))
}

mode_value <- function(x) {
  ux <- na.omit(x)
  if (length(ux) == 0L) return(NA_character_)
  names(which.max(table(ux)))
}

read_pollution_year <- function() {
  pollution_col_types <- cols(zip = col_character(), .default = col_guess())

  pm25_year <- read_csv(
    file.path(pollution_dir, "all_pm25_zip.csv.gz"),
    show_col_types = FALSE,
    col_types = pollution_col_types
  ) %>%
    filter(year >= analysis_start_year, year <= analysis_end_year) %>%
    group_by(zip = clean_zip(zip), year) %>%
    summarise(
      pm25_n_months = sum(!is.na(pm25_ug_m3)),
      pm25_ug_m3 = mean(pm25_ug_m3, na.rm = TRUE),
      pm25_source = mode_value(value_source),
      .groups = "drop"
    )

  o3_year <- read_csv(
    file.path(pollution_dir, "all_o3_zip.csv.gz"),
    show_col_types = FALSE,
    col_types = pollution_col_types
  ) %>%
    filter(year >= analysis_start_year, year <= analysis_end_year) %>%
    group_by(zip = clean_zip(zip), year) %>%
    summarise(
      o3_n_months = sum(!is.na(o3_ppb)),
      o3_ppb = mean(o3_ppb, na.rm = TRUE),
      o3_source = mode_value(value_source),
      .groups = "drop"
    )

  no2_year <- read_csv(
    file.path(pollution_dir, "all_no2_zip.csv.gz"),
    show_col_types = FALSE,
    col_types = pollution_col_types
  ) %>%
    filter(year >= analysis_start_year, year <= analysis_end_year) %>%
    transmute(zip = clean_zip(zip), year, no2 = no2, no2_source = value_source)

  pm25_year %>%
    inner_join(o3_year, by = c("zip", "year")) %>%
    inner_join(no2_year, by = c("zip", "year")) %>%
    filter(
      !is.na(pm25_ug_m3), !is.na(o3_ppb), !is.na(no2),
      pm25_n_months >= 12L, o3_n_months >= 12L
    ) %>%
    mutate(
      pm25_5ug = pm25_ug_m3 / 5,
      o3_10ppb = o3_ppb / 10,
      no2_10unit = no2 / 10
    )
}

log_msg("Reading candidate ZIP supplement")
candidate_zip <- read_sas(file.path(supp_dir, "canzip2506.sas7bdat")) %>%
  transmute(PERS_ID, PX_ID, candidate_zip = clean_zip(CAN_PERM_ZIP)) %>%
  distinct(PERS_ID, PX_ID, .keep_all = TRUE)

log_msg("Building annual pollution table")
pollution_year <- read_pollution_year()

log_msg("Reading candidate files")
candidate_data <- candidate_files %>%
  rowwise() %>%
  mutate(data = list(
    read_sas(
      path,
      col_select = any_of(c(
        "PERS_ID", "PX_ID", "WL_ORG", "CAN_LISTING_DT", "CAN_ACTIVATE_DT",
        "CAN_SOURCE", "CAN_REM_DT", "CAN_REM_CD", "CAN_DEATH_DT", "CAN_ENDWLFU",
        "CAN_AGE_AT_LISTING", "CAN_GENDER", "CAN_RACE_SRTR", "CAN_ETHNICITY_SRTR",
        "CAN_BMI", "CAN_LISTING_CTR_CD"
      ))
    ) %>%
      mutate(candidate_group = candidate_group)
  )) %>%
  ungroup() %>%
  select(data) %>%
  tidyr::unnest(data)

log_msg("Constructing baseline cohort")
cohort <- candidate_data %>%
  left_join(candidate_zip, by = c("PERS_ID", "PX_ID")) %>%
  mutate(
    waitlist_row_id = row_number(),
    index_date = coalesce(CAN_ACTIVATE_DT, CAN_LISTING_DT),
    index_year = as.integer(format(index_date, "%Y")),
    death_event = CAN_REM_CD == 8,
    transplant_event = CAN_REM_CD %in% c(4, 15, 18, 19, 14, 22, 23),
    event_date = if_else(death_event, coalesce(CAN_DEATH_DT, CAN_REM_DT), CAN_REM_DT),
    censor_date = coalesce(CAN_REM_DT, CAN_ENDWLFU),
    end_date = coalesce(event_date, censor_date),
    followup_days = as.numeric(end_date - index_date),
    event_type = case_when(
      death_event ~ 1L,
      transplant_event ~ 2L,
      TRUE ~ 0L
    ),
    age = CAN_AGE_AT_LISTING,
    sex = factor(if_else(CAN_GENDER %in% c("M", "F"), CAN_GENDER, NA_character_)),
    race_srtr = factor(na_if(CAN_RACE_SRTR, "")),
    ethnicity_srtr = factor(na_if(CAN_ETHNICITY_SRTR, "")),
    wl_org = factor(WL_ORG)
  ) %>%
  filter(
    index_year >= analysis_start_year,
    index_year <= analysis_end_year,
    CAN_SOURCE %in% c("A", "R"),
    !is.na(index_date),
    !is.na(end_date),
    !is.na(followup_days),
    followup_days >= 0
  ) %>%
  mutate(
    followup_days = pmax(followup_days, 0.5),
    death_event = as.integer(event_type == 1L),
    transplant_event = as.integer(event_type == 2L)
  ) %>%
  left_join(pollution_year, by = c("candidate_zip" = "zip", "index_year" = "year")) %>%
  filter(
    !is.na(pm25_ug_m3), !is.na(o3_ppb), !is.na(no2),
    complete.cases(age, sex, race_srtr, ethnicity_srtr, CAN_BMI, wl_org, index_year)
  ) %>%
  mutate(index_year_factor = factor(index_year))

write_csv(
  cohort %>%
    summarise(
      n_registrations = n(),
      n_people = n_distinct(PERS_ID),
      n_deaths = sum(event_type == 1L),
      n_transplants = sum(event_type == 2L),
      n_other_censored = sum(event_type == 0L),
      median_followup_days = median(followup_days)
    ),
  file.path(out_dir, "competing_risk_cohort_summary.csv")
)

log_msg("Fitting cause-specific Cox models")
cox_formula <- Surv(followup_days, death_event) ~
  pm25_5ug + o3_10ppb + no2_10unit +
  ns(age, df = 4) + sex + race_srtr + ethnicity_srtr + CAN_BMI +
  strata(wl_org) + strata(index_year_factor) + cluster(PERS_ID)

cox_death <- coxph(cox_formula, data = cohort, ties = "efron", robust = TRUE)

cox_transplant <- coxph(
  update(cox_formula, Surv(followup_days, transplant_event) ~ .),
  data = cohort,
  ties = "efron",
  robust = TRUE
)

bind_rows(
  tidy(cox_death, exponentiate = TRUE, conf.int = TRUE) %>% mutate(model = "cause_specific_waitlist_death", .before = 1),
  tidy(cox_transplant, exponentiate = TRUE, conf.int = TRUE) %>% mutate(model = "cause_specific_transplant", .before = 1)
) %>%
  filter(term %in% c("pm25_5ug", "o3_10ppb", "no2_10unit")) %>%
  transmute(
    model,
    pollutant = recode(
      term,
      pm25_5ug = "PM2.5 per 5 ug/m3",
      o3_10ppb = "O3 per 10 ppb",
      no2_10unit = "NO2 per 10 units"
    ),
    hazard_ratio = estimate,
    conf_low = conf.low,
    conf_high = conf.high,
    robust_se = std.error,
    p_value = p.value
  ) %>%
  write_csv(file.path(out_dir, "cause_specific_cox_pollutant_results.csv"))

fit_fine_gray <- function(dat, model_name) {
  fg_cov <- model.matrix(
    ~ pm25_5ug + o3_10ppb + no2_10unit + ns(age, df = 4) + sex,
    data = dat
  )[, -1, drop = FALSE]

  elapsed <- system.time({
    fit <- crr(
      ftime = dat$followup_days,
      fstatus = dat$event_type,
      cov1 = fg_cov,
      failcode = 1,
      cencode = 0
    )
  })

  tibble(
    term = names(fit$coef),
    estimate = exp(fit$coef),
    conf_low = exp(fit$coef - 1.96 * sqrt(diag(fit$var))),
    conf_high = exp(fit$coef + 1.96 * sqrt(diag(fit$var))),
    p_value = 1 - pchisq((fit$coef / sqrt(diag(fit$var)))^2, df = 1),
    elapsed_seconds = unname(elapsed[["elapsed"]])
  ) %>%
    filter(term %in% c("pm25_5ug", "o3_10ppb", "no2_10unit")) %>%
    transmute(
      model = model_name,
      pollutant = recode(
        term,
        pm25_5ug = "PM2.5 per 5 ug/m3",
        o3_10ppb = "O3 per 10 ppb",
        no2_10unit = "NO2 per 10 units"
      ),
      subdistribution_hazard_ratio = estimate,
      conf_low,
      conf_high,
      p_value,
      elapsed_seconds,
      n = nrow(dat),
      deaths = sum(dat$event_type == 1L),
      competing_transplants = sum(dat$event_type == 2L)
    )
}

log_msg("Fitting monitored organ-specific Fine-Gray models")
finegray_max_n <- as.integer(Sys.getenv("FINEGRAY_MAX_N", "10000"))
fg_groups <- cohort %>%
  group_by(WL_ORG) %>%
  group_split()
pb <- txtProgressBar(min = 0, max = length(fg_groups), style = 3)
fg_results <- vector("list", length(fg_groups))
fg_timing <- vector("list", length(fg_groups))

for (i in seq_along(fg_groups)) {
  dat <- fg_groups[[i]]
  org <- as.character(unique(dat$WL_ORG))
  if (nrow(dat) > finegray_max_n) {
    set.seed(20260527 + i)
    dat <- dat %>%
      group_by(event_type) %>%
      slice_sample(prop = min(1, finegray_max_n / nrow(fg_groups[[i]]))) %>%
      ungroup()
  }
  log_msg("Fine-Gray start: ", org, " n=", nrow(dat), " deaths=", sum(dat$event_type == 1L), " max_n=", finegray_max_n)
  start_time <- Sys.time()
  fg_results[[i]] <- tryCatch(
    fit_fine_gray(dat, paste0("fine_gray_waitlist_death_", org)),
    error = function(e) {
      tibble(
        model = paste0("fine_gray_waitlist_death_", org),
        pollutant = NA_character_,
        subdistribution_hazard_ratio = NA_real_,
        conf_low = NA_real_,
        conf_high = NA_real_,
        p_value = NA_real_,
        elapsed_seconds = as.numeric(difftime(Sys.time(), start_time, units = "secs")),
        n = nrow(dat),
        deaths = sum(dat$event_type == 1L),
        competing_transplants = sum(dat$event_type == 2L),
        error = conditionMessage(e)
      )
    }
  )
  fg_timing[[i]] <- tibble(
    WL_ORG = org,
    n = nrow(dat),
    deaths = sum(dat$event_type == 1L),
    competing_transplants = sum(dat$event_type == 2L),
    elapsed_seconds = as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  )
  log_msg("Fine-Gray done: ", org, " elapsed_seconds=", round(fg_timing[[i]]$elapsed_seconds, 1))
  setTxtProgressBar(pb, i)
}
close(pb)

bind_rows(fg_results) %>%
  write_csv(file.path(out_dir, "fine_gray_pollutant_results_by_organ.csv"))

bind_rows(fg_timing) %>%
  write_csv(file.path(out_dir, "fine_gray_timing_by_organ.csv"))

log_msg("Reading status-history files")
base_dt <- as.data.table(cohort %>%
  select(
    waitlist_row_id, PERS_ID, PX_ID, WL_ORG, index_date, end_date, event_type,
    death_event, followup_days, index_year_factor, age, sex, race_srtr,
    ethnicity_srtr, CAN_BMI, pm25_5ug, o3_10ppb, no2_10unit
  ))
setkey(base_dt, PX_ID, WL_ORG)

stathist <- stathist_files %>%
  rowwise() %>%
  mutate(data = list(
    read_sas(
      path,
      col_select = any_of(c(
        "PX_ID", "WL_ORG", "CANHX_BEGIN_DT", "CANHX_END_DT", "CANHX_STAT_CD"
      ))
    )
  )) %>%
  ungroup() %>%
  select(data) %>%
  tidyr::unnest(data) %>%
  filter(!is.na(PX_ID), !is.na(WL_ORG), !is.na(CANHX_BEGIN_DT)) %>%
  as.data.table()
setkey(stathist, PX_ID, WL_ORG)

log_msg("Constructing time-varying interval data")
tv <- stathist[base_dt, allow.cartesian = TRUE, nomatch = 0L]
tv[, hx_end := fifelse(is.na(CANHX_END_DT), end_date, CANHX_END_DT)]
tv[, interval_start_date := pmax(CANHX_BEGIN_DT, index_date)]
tv[, interval_end_date := pmin(hx_end, end_date)]
tv <- tv[
  !is.na(interval_start_date) &
    !is.na(interval_end_date) &
    interval_end_date >= interval_start_date
]
tv[, tstart := as.numeric(interval_start_date - index_date)]
tv[, tstop := as.numeric(interval_end_date - index_date)]
tv[, tstop := pmax(tstop, tstart + 0.5)]
tv[, tv_death_event := as.integer(death_event == 1L & interval_end_date >= end_date)]
tv[, current_inactive := as.integer(!is.na(CANHX_STAT_CD) & CANHX_STAT_CD %% 1000 == 999)]
tv[, current_status_unknown := as.integer(is.na(CANHX_STAT_CD) | CANHX_STAT_CD == 0)]
tv <- tv[tstop <= followup_days + 0.5]

write_csv(
  tibble(
    n_intervals = nrow(tv),
    n_registrations = uniqueN(tv$waitlist_row_id),
    n_people = uniqueN(tv$PERS_ID),
    n_death_events = sum(tv$tv_death_event),
    pct_interval_inactive = mean(tv$current_inactive),
    pct_interval_status_unknown = mean(tv$current_status_unknown)
  ),
  file.path(out_dir, "timevarying_interval_summary.csv")
)

log_msg("Fitting time-varying Cox model")
tv_fit <- coxph(
  Surv(tstart, tstop, tv_death_event) ~
    pm25_5ug + o3_10ppb + no2_10unit + current_inactive + current_status_unknown +
    ns(age, df = 4) + sex + race_srtr + ethnicity_srtr + CAN_BMI +
    strata(WL_ORG) + strata(index_year_factor) + cluster(PERS_ID),
  data = tv,
  ties = "efron",
  robust = TRUE
)

tidy(tv_fit, exponentiate = TRUE, conf.int = TRUE) %>%
  filter(term %in% c("pm25_5ug", "o3_10ppb", "no2_10unit", "current_inactive", "current_status_unknown")) %>%
  transmute(
    model = "time_varying_cox_with_current_status",
    term = recode(
      term,
      pm25_5ug = "PM2.5 per 5 ug/m3",
      o3_10ppb = "O3 per 10 ppb",
      no2_10unit = "NO2 per 10 units",
      current_inactive = "Current inactive waitlist status",
      current_status_unknown = "Current status unknown/0"
    ),
    hazard_ratio = estimate,
    conf_low = conf.low,
    conf_high = conf.high,
    robust_se = std.error,
    p_value = p.value
  ) %>%
  write_csv(file.path(out_dir, "timevarying_cox_results.csv"))

log_msg("Wrote competing-risk and time-varying outputs to ", out_dir)
