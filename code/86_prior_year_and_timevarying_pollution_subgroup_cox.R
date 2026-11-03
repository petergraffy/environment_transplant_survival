#!/usr/bin/env Rscript

runtime_source <- file.path("code", "r_runtime.R")
if (!file.exists(runtime_source)) runtime_source <- "r_runtime.R"
if (file.exists(runtime_source)) {
  source(runtime_source)
  ensure_user_library()
}

source(file.path("code", "rolling_prior_pollution.R"))

suppressPackageStartupMessages({
  library(arrow)
  library(broom)
  library(data.table)
  library(dplyr)
  library(haven)
  library(readr)
  library(scales)
  library(stringr)
  library(survival)
  library(tibble)
})

saf_paths_source <- file.path("code", "saf_paths.R")
if (!file.exists(saf_paths_source)) saf_paths_source <- "saf_paths.R"
source(saf_paths_source)
saf_paths <- get_saf_paths(release = "q1_2026")
pubsaf_dir <- saf_paths$pubsaf_dir

analysis_path <- file.path(
  "output",
  "primary_waitlist_period_pollution_cox",
  "primary_waitlist_period_pollution_analysis_dataset.csv.gz"
)
community_path <- file.path("data", "processed", "community", "zcta_acs_community_covariates_2005_2023.csv.gz")
release_dir <- file.path("data", "release")
annual_pollution_dir <- file.path(release_dir, "air_pollution_zcta_parquet")
tv_cache_dir <- file.path("output", "timevarying_pollution_cox_svi", "cache")
out_dir <- file.path("output", "prior_year_and_timevarying_pollution_subgroup_cox")
model_dir <- file.path(out_dir, "model_results")
dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tv_cache_dir, recursive = TRUE, showWarnings = FALSE)

organ_labels <- c(HR = "Heart", KI = "Kidney", LI = "Liver", LU = "Lung")
target_organs <- names(organ_labels)

log_msg <- function(...) {
  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
  flush.console()
}

clean_zip <- function(x) {
  z <- str_extract(as.character(x), "[0-9]+")
  z <- if_else(nchar(z) > 5L, substr(z, 1, 5), z)
  if_else(is.na(z), NA_character_, str_pad(z, 5, side = "left", pad = "0"))
}

format_race_group <- function(x) {
  case_when(
    x == "WHITE" ~ "White",
    x == "BLACK" ~ "Black",
    x == "ASIAN" ~ "Asian",
    TRUE ~ "Other/Unknown"
  )
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

parquet_files <- function(path) {
  list.files(path, pattern = "[.]parquet$", full.names = TRUE)
}

make_complete_acs_svi_proxy <- function(path) {
  community <- read_csv(path, show_col_types = FALSE) %>%
    mutate(zip = clean_zip(zip), analysis_year = as.integer(analysis_year))

  community_2023 <- community %>%
    filter(analysis_year == 2023L) %>%
    select(-analysis_year)

  community <- bind_rows(
    community,
    community_2023 %>% mutate(analysis_year = 2024L),
    community_2023 %>% mutate(analysis_year = 2025L)
  )

  vulnerability_vars <- c(
    "pct_poverty", "pct_unemployed", "pct_no_vehicle", "pct_nonwhite",
    "median_household_income", "pct_bachelor_plus"
  )

  community %>%
    group_by(analysis_year) %>%
    mutate(across(
      all_of(vulnerability_vars),
      ~ if_else(is.na(.x), median(.x, na.rm = TRUE), .x)
    )) %>%
    mutate(
      svi_poverty_rank = percent_rank(pct_poverty),
      svi_unemployed_rank = percent_rank(pct_unemployed),
      svi_no_vehicle_rank = percent_rank(pct_no_vehicle),
      svi_nonwhite_rank = percent_rank(pct_nonwhite),
      svi_low_income_rank = percent_rank(-median_household_income),
      svi_low_education_rank = percent_rank(-pct_bachelor_plus),
      zcta_svi_proxy = rowMeans(
        cbind(
          svi_poverty_rank,
          svi_unemployed_rank,
          svi_no_vehicle_rank,
          svi_nonwhite_rank,
          svi_low_income_rank,
          svi_low_education_rank
        ),
        na.rm = TRUE
      )
    ) %>%
    ungroup() %>%
    select(zip, analysis_year, zcta_svi_proxy)
}


daily_to_monthly_cache <- function(path, value_col, out_col, cache_file) {
  read_daily_pollution_aggregate(path, value_col, out_col, "monthly", cache_file)
}


read_exposure_tables <- function() {
  pm25_monthly <- daily_to_monthly_cache(
    file.path(release_dir, "lghap_pm25_zcta_daily_parquet"),
    "pm25_ug_m3",
    "pm25_interval_ug_m3",
    file.path(tv_cache_dir, "pm25_daily_derived_monthly_zcta.csv.gz")
  )
  o3_monthly <- daily_to_monthly_cache(
    file.path(release_dir, "o3_zcta_daily_parquet"),
    "o3_ppb",
    "o3_interval_ppb",
    file.path(tv_cache_dir, "o3_daily_derived_monthly_zcta.csv.gz")
  )
  no2_monthly <- open_dataset(parquet_files(file.path(release_dir, "no2_zcta_monthly_parquet"))) %>%
    transmute(zip = zip, year = year, month = month, no2_interval_ppb = no2_ppbv) %>%
    collect() %>%
    mutate(zip = clean_zip(zip), year = as.integer(year), month = as.integer(month))
  no2_annual <- read_parquet(file.path(annual_pollution_dir, "air_pollution_zcta_no2_annual_2005_2025.parquet")) %>%
    transmute(zip = clean_zip(zip), year = as.integer(year), no2_annual_ppb = no2)

  list(
    pm25 = as.data.table(pm25_monthly),
    o3 = as.data.table(o3_monthly),
    no2_monthly = as.data.table(no2_monthly),
    no2_annual = as.data.table(no2_annual)
  )
}

make_month_intervals <- function(cohort_dt, exposure_end_date) {
  intervals <- trim_monthly_cohort(cohort_dt)
  intervals[, analysis_end_date := pmin(observed_end_date, exposure_end_date)]
  intervals <- intervals[analysis_end_date >= index_date]
  out <- expand_monthly_cohort(intervals)
  out[, interval_start_date := pmax(index_date, month_start)]
  out[, interval_end_date := pmin(analysis_end_date, month_end)]
  out <- out[interval_end_date >= interval_start_date]
  out[, `:=`(
    year = as.integer(format(month_start, "%Y")),
    month = as.integer(format(month_start, "%m")),
    tstart = as.numeric(interval_start_date - index_date),
    tstop = as.numeric(interval_end_date - index_date) + 1,
    tv_adverse_event = as.integer(adverse_event == 1L & observed_end_date <= exposure_end_date & interval_end_date >= observed_end_date)
  )]
  out[, tstop := pmax(tstop, tstart + 0.5)]
  out[, kidney_dialysis_years_tv := kidney_dialysis_years]
  out[WL_ORG == "KI" & kidney_no_dialysis_time == 0L & !is.na(kidney_dialysis_years), kidney_dialysis_years_tv := kidney_dialysis_years + tstart / 365.25]
  out[, organ_score_tv := organ_score]
  out[, interval_row_id := .I]
  out
}

make_time_updated_score_sources <- function() {
  log_msg("Reading time-updated liver MELD and heart US-CRS inputs")
  liver_updates <- read_sas(
    file.path(pubsaf_dir, "stathist_liin.sas7bdat"),
    col_select = any_of(c("PX_ID", "WL_ORG", "CANHX_BEGIN_DT", "CANHX_SRTR_LAB_MELD", "CANHX_OPTN_LAB_MELD"))
  ) %>%
    filter(WL_ORG == "LI", !is.na(PX_ID), !is.na(CANHX_BEGIN_DT)) %>%
    transmute(
      PX_ID,
      score_date = as.Date(CANHX_BEGIN_DT),
      liver_score_update = coalesce(CANHX_SRTR_LAB_MELD, CANHX_OPTN_LAB_MELD),
      liver_score_update = if_else(!is.na(liver_score_update) & liver_score_update > 1000, liver_score_update - 6200, liver_score_update)
    ) %>%
    filter(is.finite(liver_score_update)) %>%
    as.data.table()

  heart_updates <- read_sas(
    file.path(pubsaf_dir, "statjust_hr1a.sas7bdat"),
    col_select = any_of(c(
      "PX_ID", "WL_ORG", "CANHX_CHG_DT", "CANHX_RVAD_TYPE", "CANHX_ECMO",
      "CANHX_LVAD_TYPE", "CANHX_LAB_SERUM_CREAT", "CANHX_LAB_BILI",
      "CANHX_LAB_ALBUMIN", "CANHX_LAB_SODIUM", "CANHX_LAB_BNP"
    ))
  ) %>%
    filter(WL_ORG == "HR", !is.na(PX_ID), !is.na(CANHX_CHG_DT)) %>%
    transmute(
      PX_ID,
      score_date = as.Date(CANHX_CHG_DT),
      hr_short_mcs_update = as.integer(flag_yes(CANHX_ECMO) == 1L | flag_yes(CANHX_RVAD_TYPE) == 1L),
      hr_durable_lvad_update = flag_yes(CANHX_LVAD_TYPE),
      hr_creatinine_update = as.numeric(CANHX_LAB_SERUM_CREAT),
      hr_bilirubin_update = as.numeric(CANHX_LAB_BILI),
      hr_albumin_update = as.numeric(CANHX_LAB_ALBUMIN),
      hr_sodium_update = as.numeric(CANHX_LAB_SODIUM),
      hr_bnp_update = as.numeric(CANHX_LAB_BNP)
    ) %>%
    as.data.table()

  if (nrow(heart_updates) > 0L) {
    heart_updates[, `:=`(
      hr_creatinine_update = median_impute_num(hr_creatinine_update),
      hr_bilirubin_update = median_impute_num(hr_bilirubin_update),
      hr_albumin_update = median_impute_num(hr_albumin_update),
      hr_sodium_update = median_impute_num(hr_sodium_update),
      hr_bnp_update = median_impute_num(hr_bnp_update)
    )]
  }

  setkey(liver_updates, PX_ID, score_date)
  setkey(heart_updates, PX_ID, score_date)
  list(liver = liver_updates, heart = heart_updates)
}

add_time_updated_scores <- function(intervals, score_sources) {
  intervals[, age_interval := age + tstart / 365.25]

  li_intervals <- intervals[WL_ORG == "LI", .(interval_row_id, PX_ID, interval_start_date)]
  if (nrow(li_intervals) > 0L && nrow(score_sources$liver) > 0L) {
    setkey(li_intervals, PX_ID, interval_start_date)
    li_join <- score_sources$liver[li_intervals, on = c("PX_ID", "score_date" = "interval_start_date"), roll = Inf]
    intervals[li_join, organ_score_tv := i.liver_score_update, on = "interval_row_id"]
  }

  hr_intervals <- intervals[WL_ORG == "HR", .(interval_row_id, PX_ID, interval_start_date, age_interval, sex)]
  if (nrow(hr_intervals) > 0L && nrow(score_sources$heart) > 0L) {
    setkey(hr_intervals, PX_ID, interval_start_date)
    hr_join <- score_sources$heart[hr_intervals, on = c("PX_ID", "score_date" = "interval_start_date"), roll = Inf]
    hr_join[, hr_egfr_update := calc_egfr_2021(hr_creatinine_update, age_interval, sex)]
    hr_join[, heart_score_update := 1.02 * hr_short_mcs_update +
      0.55 * safe_log1p(hr_bilirubin_update) -
      0.01 * hr_egfr_update +
      0.40 * safe_log1p(hr_bnp_update) -
      0.63 * hr_albumin_update -
      0.07 * hr_sodium_update -
      1.12 * hr_durable_lvad_update]
    intervals[hr_join, organ_score_tv := i.heart_score_update, on = "interval_row_id"]
  }

  intervals[is.na(organ_score_tv), organ_score_tv := organ_score]
  intervals
}

subgroup_specs <- list(
  age_group = list(label = "Age group", baseline_adjust = c("sex", "race_group", "zcta_svi_proxy"), tv_adjust = c("sex", "race_group", "zcta_svi_proxy")),
  sex = list(label = "Sex", baseline_adjust = c("age", "race_group", "zcta_svi_proxy"), tv_adjust = c("age_interval", "race_group", "zcta_svi_proxy")),
  race_group = list(label = "Race", baseline_adjust = c("age", "sex", "zcta_svi_proxy"), tv_adjust = c("age_interval", "sex", "zcta_svi_proxy"))
)

subgroup_levels <- list(
  age_group = c("<18", "18-39", "40-59", "60+"),
  sex = c("Female", "Male"),
  race_group = c("White", "Black", "Asian")
)

tv_organ_terms <- function(org) {
  if (org == "KI") {
    c("kidney_no_dialysis_time", "kidney_dialysis_years_tv", "kidney_diabetes")
  } else {
    "organ_score_tv"
  }
}

fit_baseline_subgroup <- function(dat, org, pollutant, exposure_term, exposure_label, subgroup_var, subgroup_level, subgroup_spec) {
  adjust_terms <- subgroup_spec$baseline_adjust
  vars_needed <- c("followup_days", "adverse_event", exposure_term, subgroup_var, adjust_terms, "listing_center")
  model_dat <- dat %>%
    filter(
      WL_ORG == org,
      .data[[subgroup_var]] == subgroup_level,
      complete.cases(across(all_of(unique(vars_needed)))),
      followup_days > 0
    ) %>%
    mutate(.followup_days = pmax(followup_days, 0.5)) %>%
    droplevels()

  source_n <- nrow(model_dat)
  source_events <- sum(model_dat$adverse_event, na.rm = TRUE)
  if (source_n == 0L || source_events < 50L || n_distinct(model_dat$listing_center) < 2L) {
    return(tibble(
      model = "baseline_prior_year",
      organ = org,
      organ_label = recode(org, !!!organ_labels),
      pollutant = pollutant,
      exposure_label = exposure_label,
      subgroup = subgroup_var,
      subgroup_label = subgroup_spec$label,
      subgroup_level = subgroup_level,
      n = source_n,
      people = n_distinct(model_dat$PERS_ID),
      centers = n_distinct(model_dat$listing_center),
      adverse_events = source_events,
      hazard_ratio = NA_real_,
      conf_low = NA_real_,
      conf_high = NA_real_,
      p_value = NA_real_,
      adjustment_set = paste(c(adjust_terms, "strata(listing_center)"), collapse = " + "),
      error = "Insufficient sample/events/centers"
    ))
  }

  log_msg("Baseline subgroup Cox ", org, " | ", pollutant, " | ", subgroup_var, "=", subgroup_level, " n=", source_n, " adverse=", source_events)
  form <- as.formula(paste(
    "Surv(.followup_days, adverse_event) ~",
    paste(c(exposure_term, adjust_terms, "strata(listing_center)"), collapse = " + ")
  ))
  fit <- coxph(form, data = model_dat, ties = "efron", x = FALSE, y = FALSE)

  tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(term == exposure_term) %>%
    transmute(
      model = "baseline_prior_year",
      organ = org,
      organ_label = recode(org, !!!organ_labels),
      pollutant = pollutant,
      exposure_label = exposure_label,
      subgroup = subgroup_var,
      subgroup_label = subgroup_spec$label,
      subgroup_level = subgroup_level,
      n = source_n,
      people = n_distinct(model_dat$PERS_ID),
      centers = n_distinct(model_dat$listing_center),
      adverse_events = source_events,
      hazard_ratio = estimate,
      conf_low = conf.low,
      conf_high = conf.high,
      p_value = p.value,
      adjustment_set = paste(c(adjust_terms, "strata(listing_center)"), collapse = " + "),
      error = NA_character_
    )
}

fit_tv_subgroup <- function(intervals, org, pollutant, exposure_term, exposure_label, subgroup_var, subgroup_level, subgroup_spec) {
  adjust_terms <- c(subgroup_spec$tv_adjust, tv_organ_terms(org))
  vars_needed <- c("tstart", "tstop", "tv_adverse_event", exposure_term, subgroup_var, adjust_terms, "listing_center")
  model_dat <- intervals[
    WL_ORG == org &
      get(subgroup_var) == subgroup_level &
      complete.cases(intervals[, ..vars_needed]) &
      tstop > tstart
  ]

  source_n <- nrow(model_dat)
  source_events <- sum(model_dat$tv_adverse_event, na.rm = TRUE)
  if (source_n == 0L || source_events < 50L || uniqueN(model_dat$listing_center) < 2L) {
    return(tibble(
      model = "time_varying",
      organ = org,
      organ_label = recode(org, !!!organ_labels),
      pollutant = pollutant,
      exposure_label = exposure_label,
      subgroup = subgroup_var,
      subgroup_label = subgroup_spec$label,
      subgroup_level = subgroup_level,
      n = source_n,
      people = uniqueN(model_dat$PERS_ID),
      centers = uniqueN(model_dat$listing_center),
      adverse_events = source_events,
      hazard_ratio = NA_real_,
      conf_low = NA_real_,
      conf_high = NA_real_,
      p_value = NA_real_,
      adjustment_set = paste(c(adjust_terms, "strata(listing_center)"), collapse = " + "),
      error = "Insufficient sample/events/centers"
    ))
  }

  log_msg("Time-varying subgroup Cox ", org, " | ", pollutant, " | ", subgroup_var, "=", subgroup_level, " intervals=", source_n, " adverse=", source_events)
  form <- as.formula(paste(
    "Surv(tstart, tstop, tv_adverse_event) ~",
    paste(c(exposure_term, adjust_terms, "strata(listing_center)"), collapse = " + ")
  ))
  fit <- coxph(form, data = model_dat, ties = "efron", x = FALSE, y = FALSE)

  tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(term == exposure_term) %>%
    transmute(
      model = "time_varying",
      organ = org,
      organ_label = recode(org, !!!organ_labels),
      pollutant = pollutant,
      exposure_label = exposure_label,
      subgroup = subgroup_var,
      subgroup_label = subgroup_spec$label,
      subgroup_level = subgroup_level,
      n = source_n,
      people = uniqueN(model_dat$PERS_ID),
      centers = uniqueN(model_dat$listing_center),
      adverse_events = source_events,
      hazard_ratio = estimate,
      conf_low = conf.low,
      conf_high = conf.high,
      p_value = p.value,
      adjustment_set = paste(c(adjust_terms, "strata(listing_center)"), collapse = " + "),
      error = NA_character_
    )
}

log_msg("Reading primary deduplicated cohort")
analysis_dat <- read_csv(analysis_path, show_col_types = FALSE) %>%
  mutate(
    candidate_zip = clean_zip(candidate_zip),
    index_date = as.Date(index_date),
    observed_end_date = as.Date(observed_end_date),
    listing_year_int = as.integer(as.character(listing_year)),
    followup_days = as.numeric(observed_end_date - index_date),
    sex = recode(as.character(sex), F = "Female", M = "Male", .default = as.character(sex)),
    sex = factor(sex, levels = c("Female", "Male")),
    race_group = factor(format_race_group(race), levels = c("White", "Black", "Asian", "Other/Unknown")),
    race = factor(race),
    age_group = cut(age, breaks = c(-Inf, 17, 39, 59, Inf), labels = c("<18", "18-39", "40-59", "60+"), right = TRUE),
    age_group = factor(age_group, levels = c("<18", "18-39", "40-59", "60+")),
    listing_center = factor(listing_center)
  ) %>%
  filter(WL_ORG %in% target_organs, observed_end_date >= index_date, !is.na(candidate_zip))

log_msg("Attaching ACS-derived ZCTA SVI proxy")
svi <- make_complete_acs_svi_proxy(community_path)
analysis_dat <- analysis_dat %>%
  left_join(svi, by = c("candidate_zip" = "zip", "listing_year_int" = "analysis_year"), suffix = c("", "_community")) %>%
  mutate(zcta_svi_proxy = coalesce(zcta_svi_proxy, zcta_svi_proxy_community)) %>%
  select(-any_of("zcta_svi_proxy_community"))

log_msg("Attaching prior-year pollution")
prior_pollution <- read_rolling_prior_pollution(analysis_dat)
analysis_dat <- analysis_dat %>%
  left_join(prior_pollution$pm25, by = c("candidate_zip" = "zip", "index_date" = "index_date")) %>%
  left_join(prior_pollution$o3, by = c("candidate_zip" = "zip", "index_date" = "index_date")) %>%
  left_join(prior_pollution$no2, by = c("candidate_zip" = "zip", "index_date" = "index_date")) %>%
  mutate(
    pm25_prior_5ug = pm25_prior_ug_m3 / 5,
    o3_prior_10ppb = o3_prior_ppb / 10,
    no2_prior_10ppb = no2_prior_ppb / 10
  )

pollutant_specs <- tribble(
  ~pollutant, ~baseline_term, ~tv_term, ~label, ~tv_label, ~exposure_end_date,
  "pm25", "pm25_prior_5ug", "pm25_interval_5ug", "PM2.5 per 5 ug/m3", "Monthly PM2.5 per 5 ug/m3", daily_pollution_end_date,
  "o3", "o3_prior_10ppb", "o3_interval_10ppb", "O3 per 10 ppb", "Monthly O3 per 10 ppb", daily_pollution_end_date,
  "no2", "no2_prior_10ppb", "no2_interval_10ppb", "NO2 per 10 ppb", "Monthly/annual NO2 per 10 ppb", as.Date("2025-12-31")
)

baseline_results <- list()
for (org in target_organs) {
  for (i in seq_len(nrow(pollutant_specs))) {
    spec <- pollutant_specs[i, ]
    for (subgroup_var in names(subgroup_specs)) {
      subgroup_spec <- subgroup_specs[[subgroup_var]]
      for (level in subgroup_levels[[subgroup_var]]) {
        baseline_results[[length(baseline_results) + 1L]] <- fit_baseline_subgroup(
          analysis_dat, org, spec$pollutant, spec$baseline_term, spec$label, subgroup_var, level, subgroup_spec
        )
      }
    }
  }
}

baseline_results <- bind_rows(baseline_results)
write_csv(baseline_results, file.path(out_dir, "prior_year_pollution_subgroup_cox_results.csv"))

log_msg("Reading time-varying exposure tables")
exposure_tables <- read_exposure_tables()
for (nm in names(exposure_tables)) {
  if (nm %in% c("pm25", "o3", "no2_monthly")) setkey(exposure_tables[[nm]], zip, year, month)
  if (nm == "no2_annual") setkey(exposure_tables[[nm]], zip, year)
}

score_sources <- make_time_updated_score_sources()
analysis_dt <- as.data.table(analysis_dat)
tv_results <- list()
tv_summaries <- list()

for (i in seq_len(nrow(pollutant_specs))) {
  spec <- pollutant_specs[i, ]
  log_msg("Building monthly time-varying subgroup intervals for ", spec$pollutant)
  intervals <- make_month_intervals(analysis_dt, spec$exposure_end_date)
  intervals <- add_time_updated_scores(intervals, score_sources)
  intervals[, pollutant := spec$pollutant]

  if (spec$pollutant == "pm25") {
    intervals <- exposure_tables$pm25[intervals, on = c("zip" = "candidate_zip", "year", "month")]
    intervals[, pm25_interval_5ug := pm25_interval_ug_m3 / 5]
  } else if (spec$pollutant == "o3") {
    intervals <- exposure_tables$o3[intervals, on = c("zip" = "candidate_zip", "year", "month")]
    intervals[, o3_interval_10ppb := o3_interval_ppb / 10]
  } else {
    intervals <- exposure_tables$no2_monthly[intervals, on = c("zip" = "candidate_zip", "year", "month")]
    intervals <- exposure_tables$no2_annual[intervals, on = c("zip", "year")]
    intervals[, no2_interval_ppb_final := fifelse(!is.na(no2_interval_ppb), no2_interval_ppb, no2_annual_ppb)]
    intervals[, no2_interval_10ppb := no2_interval_ppb_final / 10]
  }

  tv_summaries[[spec$pollutant]] <- intervals[
    ,
    .(
      intervals = .N,
      candidate_organ_episodes = uniqueN(waitlist_row_id),
      people = uniqueN(PERS_ID),
      adverse_events = sum(tv_adverse_event),
      complete_exposure_intervals = sum(!is.na(get(spec$tv_term))),
      min_year = min(year, na.rm = TRUE),
      max_year = max(year, na.rm = TRUE)
    ),
    by = .(pollutant, WL_ORG)
  ]

  for (org in target_organs) {
    for (subgroup_var in names(subgroup_specs)) {
      subgroup_spec <- subgroup_specs[[subgroup_var]]
      for (level in subgroup_levels[[subgroup_var]]) {
        tv_results[[length(tv_results) + 1L]] <- fit_tv_subgroup(
          intervals, org, spec$pollutant, spec$tv_term, spec$tv_label, subgroup_var, level, subgroup_spec
        )
      }
    }
  }

  rm(intervals)
  invisible(gc())
}

tv_results <- bind_rows(tv_results)
write_csv(tv_results, file.path(out_dir, "timevarying_pollution_subgroup_cox_results.csv"))
write_csv(bind_rows(tv_summaries), file.path(out_dir, "timevarying_pollution_subgroup_interval_summary.csv"))

combined <- bind_rows(baseline_results, tv_results) %>%
  mutate(
    hazard_ratio_ci = if_else(
      is.na(hazard_ratio),
      NA_character_,
      sprintf("%.2f (%.2f-%.2f)", hazard_ratio, conf_low, conf_high)
    ),
    p_value_display = case_when(
      is.na(p_value) ~ NA_character_,
      p_value < 0.001 ~ "<.001",
      TRUE ~ sprintf("%.3f", p_value)
    )
  )
write_csv(combined, file.path(out_dir, "prior_year_and_timevarying_pollution_subgroup_cox_results.csv"))

log_msg("Wrote prior-year and time-varying subgroup Cox outputs to ", normalizePath(out_dir, winslash = "/"))
