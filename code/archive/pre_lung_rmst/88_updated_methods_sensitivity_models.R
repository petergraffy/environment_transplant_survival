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
  library(tidyr)
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
out_dir <- file.path("output", "updated_methods_sensitivity_models")
model_dir <- file.path(out_dir, "model_results")
table_dir <- file.path("output", "tables", "updated_methods_sensitivity_models")
dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
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

format_p <- function(p) {
  case_when(
    is.na(p) ~ "",
    p < 0.001 ~ "<.001",
    TRUE ~ sprintf("%.3f", p)
  )
}

format_est_ci <- function(est, low, high, digits = 2) {
  sprintf(paste0("%.", digits, "f (%.", digits, "f-%.", digits, "f)"), est, low, high)
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

attach_monthly_exposures <- function(intervals, exposure_tables) {
  pm25 <- copy(exposure_tables$pm25)
  o3 <- copy(exposure_tables$o3)
  no2_monthly <- copy(exposure_tables$no2_monthly)
  no2_annual <- copy(exposure_tables$no2_annual)
  setnames(pm25, "zip", "candidate_zip")
  setnames(o3, "zip", "candidate_zip")
  setnames(no2_monthly, "zip", "candidate_zip")
  setnames(no2_annual, "zip", "candidate_zip")

  intervals <- pm25[intervals, on = c("candidate_zip", "year", "month")]
  intervals[, pm25_interval_5ug := pm25_interval_ug_m3 / 5]
  intervals <- o3[intervals, on = c("candidate_zip", "year", "month")]
  intervals[, o3_interval_10ppb := o3_interval_ppb / 10]
  intervals <- no2_monthly[intervals, on = c("candidate_zip", "year", "month")]
  intervals <- no2_annual[intervals, on = c("candidate_zip", "year")]
  intervals[, no2_interval_ppb_final := fifelse(!is.na(no2_interval_ppb), no2_interval_ppb, no2_annual_ppb)]
  intervals[, no2_interval_10ppb := no2_interval_ppb_final / 10]
  intervals
}

tv_adjustment <- function(org, include_center = TRUE) {
  terms <- c("age", "sex", "race", "zcta_svi_proxy")
  vars <- terms
  if (org == "KI") {
    terms <- c(terms, "kidney_no_dialysis_time", "kidney_dialysis_years_tv", "kidney_diabetes")
    vars <- c(vars, "kidney_no_dialysis_time", "kidney_dialysis_years_tv", "kidney_diabetes")
  } else {
    terms <- c(terms, "organ_score_tv")
    vars <- c(vars, "organ_score_tv")
  }
  if (include_center) {
    terms <- c(terms, "strata(listing_center)")
    vars <- c(vars, "listing_center")
  }
  list(terms = terms, vars = vars)
}

fit_cox_terms <- function(fit, exposure_terms) {
  tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(term %in% exposure_terms) %>%
    transmute(
      term = term,
      hazard_ratio = estimate,
      conf_low = conf.low,
      conf_high = conf.high,
      p_value = p.value
    )
}

fit_baseline_model <- function(dat, org, model_name, exposure_terms, exposure_labels,
                               include_center, include_multi_organ = FALSE,
                               exclude_multi_organ = FALSE) {
  adjustment_terms <- c("age", "sex", "race", "zcta_svi_proxy")
  if (include_multi_organ) adjustment_terms <- c(adjustment_terms, "multi_organ_candidate")
  if (include_center) adjustment_terms <- c(adjustment_terms, "strata(listing_center)")
  vars_needed <- c("followup_days", "adverse_event", exposure_terms, "age", "sex", "race", "zcta_svi_proxy")
  if (include_multi_organ) vars_needed <- c(vars_needed, "multi_organ_candidate")
  if (include_center) vars_needed <- c(vars_needed, "listing_center")
  model_dat <- dat %>%
    filter(
      WL_ORG == org,
      !exclude_multi_organ | multi_organ_candidate == "Single-organ candidate",
      complete.cases(across(all_of(vars_needed))),
      followup_days > 0
    ) %>%
    mutate(.followup_days = pmax(followup_days, 0.5)) %>%
    droplevels()

  if (nrow(model_dat) == 0L || sum(model_dat$adverse_event) < 50L) return(tibble())

  log_msg("Baseline sensitivity Cox ", org, " | ", model_name, " n=", nrow(model_dat), " adverse=", sum(model_dat$adverse_event))
  form <- as.formula(paste(
    "Surv(.followup_days, adverse_event) ~",
    paste(c(exposure_terms, adjustment_terms), collapse = " + ")
  ))
  fit <- coxph(form, data = model_dat, ties = "efron", x = FALSE, y = FALSE)
  fit_cox_terms(fit, exposure_terms) %>%
    mutate(
      analysis = "Baseline cause-specific Cox",
      model = model_name,
      organ = org,
      organ_label = recode(org, !!!organ_labels),
      exposure = term,
      exposure_label = unname(exposure_labels[term]),
      n = nrow(model_dat),
      people = n_distinct(model_dat$PERS_ID),
      intervals = NA_integer_,
      centers = n_distinct(model_dat$listing_center),
      adverse_events = sum(model_dat$adverse_event),
      adjustment_set = paste(adjustment_terms, collapse = " + "),
      .before = term
    )
}

fit_tv_model <- function(intervals, org, model_name, exposure_terms, exposure_labels,
                         include_center, include_multi_organ = FALSE,
                         exclude_multi_organ = FALSE) {
  adj <- tv_adjustment(org, include_center = include_center)
  if (include_multi_organ) {
    adj$terms <- c("multi_organ_candidate", adj$terms)
    adj$vars <- c("multi_organ_candidate", adj$vars)
  }
  vars_needed <- c("tstart", "tstop", "tv_adverse_event", exposure_terms, adj$vars)
  model_dat <- intervals[
    WL_ORG == org &
      (!exclude_multi_organ | multi_organ_candidate == "Single-organ candidate") &
      complete.cases(intervals[, ..vars_needed]) &
      tstop > tstart
  ]

  if (nrow(model_dat) == 0L || sum(model_dat$tv_adverse_event) < 50L) return(tibble())

  log_msg("Time-varying sensitivity Cox ", org, " | ", model_name, " intervals=", nrow(model_dat), " adverse=", sum(model_dat$tv_adverse_event))
  form <- as.formula(paste(
    "Surv(tstart, tstop, tv_adverse_event) ~",
    paste(c(exposure_terms, adj$terms), collapse = " + ")
  ))
  fit <- coxph(form, data = model_dat, ties = "efron", x = FALSE, y = FALSE)
  fit_cox_terms(fit, exposure_terms) %>%
    mutate(
      analysis = "Time-varying cause-specific Cox",
      model = model_name,
      organ = org,
      organ_label = recode(org, !!!organ_labels),
      exposure = term,
      exposure_label = unname(exposure_labels[term]),
      n = uniqueN(model_dat$waitlist_row_id),
      people = uniqueN(model_dat$PERS_ID),
      intervals = nrow(model_dat),
      centers = if ("listing_center" %in% names(model_dat)) uniqueN(model_dat$listing_center) else NA_integer_,
      adverse_events = sum(model_dat$tv_adverse_event),
      adjustment_set = paste(adj$terms, collapse = " + "),
      .before = term
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
    sex = factor(sex),
    race = factor(race),
    listing_center = factor(listing_center)
  ) %>%
  filter(WL_ORG %in% target_organs, observed_end_date >= index_date, !is.na(candidate_zip))

multi_organ_people <- analysis_dat %>%
  distinct(PERS_ID, WL_ORG) %>%
  count(PERS_ID, name = "n_organ_groups") %>%
  mutate(
    multi_organ_candidate = factor(
      if_else(n_organ_groups > 1L, "Multi-organ candidate", "Single-organ candidate"),
      levels = c("Single-organ candidate", "Multi-organ candidate")
    )
  )

analysis_dat <- analysis_dat %>%
  left_join(multi_organ_people, by = "PERS_ID")

multi_organ_summary <- analysis_dat %>%
  distinct(PERS_ID, multi_organ_candidate, n_organ_groups) %>%
  count(multi_organ_candidate, n_organ_groups, name = "people") %>%
  arrange(n_organ_groups, multi_organ_candidate)
write_csv(multi_organ_summary, file.path(out_dir, "multi_organ_candidate_summary.csv"))

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

baseline_label_terms <- c(
  pm25 = "pm25_prior_5ug",
  o3 = "o3_prior_10ppb",
  no2 = "no2_prior_10ppb"
)
baseline_term_labels <- c(
  pm25_prior_5ug = "PM2.5 per 5 ug/m3",
  o3_prior_10ppb = "O3 per 10 ppb",
  no2_prior_10ppb = "NO2 per 10 ppb"
)
tv_label_terms <- c(
  pm25 = "pm25_interval_5ug",
  o3 = "o3_interval_10ppb",
  no2 = "no2_interval_10ppb"
)
tv_term_labels <- c(
  pm25_interval_5ug = "PM2.5 per 5 ug/m3",
  o3_interval_10ppb = "O3 per 10 ppb",
  no2_interval_10ppb = "NO2 per 10 ppb"
)

baseline_results <- list()
for (org in target_organs) {
  for (pollutant in names(baseline_label_terms)) {
    term <- baseline_label_terms[[pollutant]]
    baseline_results[[length(baseline_results) + 1L]] <- fit_baseline_model(
      analysis_dat,
      org,
      "Single-pollutant model restricted to single-organ candidates",
      term,
      setNames(baseline_term_labels[[term]], term),
      include_center = TRUE,
      exclude_multi_organ = TRUE
    )
    baseline_results[[length(baseline_results) + 1L]] <- fit_baseline_model(
      analysis_dat,
      org,
      "Single-pollutant model additionally adjusted for multi-organ candidate status",
      term,
      setNames(baseline_term_labels[[term]], term),
      include_center = TRUE,
      include_multi_organ = TRUE
    )
    baseline_results[[length(baseline_results) + 1L]] <- fit_baseline_model(
      analysis_dat,
      org,
      "Single-pollutant model without listing center",
      term,
      setNames(baseline_term_labels[[term]], term),
      include_center = FALSE
    )
  }
  baseline_results[[length(baseline_results) + 1L]] <- fit_baseline_model(
    analysis_dat,
    org,
    "Multipollutant model: PM2.5 + NO2",
    unname(baseline_label_terms[c("pm25", "no2")]),
    baseline_term_labels[unname(baseline_label_terms[c("pm25", "no2")])],
    include_center = TRUE
  )
  baseline_results[[length(baseline_results) + 1L]] <- fit_baseline_model(
    analysis_dat,
    org,
    "Multipollutant model: PM2.5 + NO2 + O3",
    unname(baseline_label_terms[c("pm25", "no2", "o3")]),
    baseline_term_labels[unname(baseline_label_terms[c("pm25", "no2", "o3")])],
    include_center = TRUE
  )
}
baseline_results <- bind_rows(baseline_results)
write_csv(baseline_results, file.path(out_dir, "baseline_prior_year_sensitivity_results.csv"))

log_msg("Reading time-varying exposure tables")
exposure_tables <- read_exposure_tables()
for (nm in names(exposure_tables)) {
  if (nm %in% c("pm25", "o3", "no2_monthly")) setkey(exposure_tables[[nm]], zip, year, month)
  if (nm == "no2_annual") setkey(exposure_tables[[nm]], zip, year)
}
score_sources <- make_time_updated_score_sources()
analysis_dt <- as.data.table(analysis_dat)

tv_results <- list()
tv_end_dates <- c(pm25 = daily_pollution_end_date, o3 = daily_pollution_end_date, no2 = as.Date("2025-12-31"))
for (pollutant in names(tv_label_terms)) {
  log_msg("Building time-varying intervals for no-center ", pollutant)
  intervals <- make_month_intervals(analysis_dt, tv_end_dates[[pollutant]])
  intervals <- add_time_updated_scores(intervals, score_sources)
  intervals <- attach_monthly_exposures(intervals, exposure_tables)
    for (org in target_organs) {
    tv_results[[length(tv_results) + 1L]] <- fit_tv_model(
      intervals,
      org,
      "Single-pollutant model restricted to single-organ candidates",
      tv_label_terms[[pollutant]],
      setNames(tv_term_labels[[tv_label_terms[[pollutant]]]], tv_label_terms[[pollutant]]),
      include_center = TRUE,
      exclude_multi_organ = TRUE
    )
    tv_results[[length(tv_results) + 1L]] <- fit_tv_model(
      intervals,
      org,
      "Single-pollutant model additionally adjusted for multi-organ candidate status",
      tv_label_terms[[pollutant]],
      setNames(tv_term_labels[[tv_label_terms[[pollutant]]]], tv_label_terms[[pollutant]]),
      include_center = TRUE,
      include_multi_organ = TRUE
    )
    tv_results[[length(tv_results) + 1L]] <- fit_tv_model(
      intervals,
      org,
      "Single-pollutant model without listing center",
      tv_label_terms[[pollutant]],
      setNames(tv_term_labels[[tv_label_terms[[pollutant]]]], tv_label_terms[[pollutant]]),
      include_center = FALSE
    )
  }
  rm(intervals)
  invisible(gc())
}

log_msg("Building common PM2.5-limited intervals for time-varying multipollutant models")
multi_intervals <- make_month_intervals(analysis_dt, daily_pollution_end_date)
multi_intervals <- add_time_updated_scores(multi_intervals, score_sources)
multi_intervals <- attach_monthly_exposures(multi_intervals, exposure_tables)
for (org in target_organs) {
  tv_results[[length(tv_results) + 1L]] <- fit_tv_model(
    multi_intervals,
    org,
    "Multipollutant model: PM2.5 + NO2",
    unname(tv_label_terms[c("pm25", "no2")]),
    tv_term_labels[unname(tv_label_terms[c("pm25", "no2")])],
    include_center = TRUE
  )
  tv_results[[length(tv_results) + 1L]] <- fit_tv_model(
    multi_intervals,
    org,
    "Multipollutant model: PM2.5 + NO2 + O3",
    unname(tv_label_terms[c("pm25", "no2", "o3")]),
    tv_term_labels[unname(tv_label_terms[c("pm25", "no2", "o3")])],
    include_center = TRUE
  )
}
rm(multi_intervals)
invisible(gc())

tv_results <- bind_rows(tv_results)
write_csv(tv_results, file.path(out_dir, "timevarying_sensitivity_results.csv"))

sensitivity_results <- bind_rows(baseline_results, tv_results) %>%
  select(
    analysis, model, organ, organ_label, exposure, exposure_label, n, people, intervals,
    centers, adverse_events, hazard_ratio, conf_low, conf_high, p_value,
    adjustment_set
  )
write_csv(sensitivity_results, file.path(out_dir, "updated_methods_sensitivity_results.csv"))

cox_table <- sensitivity_results %>%
  mutate(
    section = "Cox model sensitivity analyses",
    result = format_est_ci(hazard_ratio, conf_low, conf_high, digits = 2),
    p_value_display = format_p(p_value)
  ) %>%
  transmute(
    Section = section,
    Analysis = analysis,
    Model = model,
    Organ = organ_label,
    Pollutant = exposure_label,
    `N candidates` = n,
    `Adverse events` = adverse_events,
    `HR (95% CI)` = result,
    `P value` = p_value_display,
    `Adjustment set` = adjustment_set
  )

baseline_severity <- read_csv(
  file.path("output", "prior_year_pollution_baseline_severity_associations", "prior_year_pollution_baseline_severity_associations.csv"),
  show_col_types = FALSE
) %>%
  filter(measure == "mean_difference") %>%
  mutate(
    analysis = "Baseline organ severity association",
    model = "Pollution exposure as predictor of baseline organ severity",
    n_display = n,
    result = format_est_ci(estimate, conf_low, conf_high, digits = 2),
    p_value_display = format_p(p_value)
  ) %>%
  transmute(
    Section = "Pollution-organ severity associations",
    Analysis = analysis,
    Model = model,
    Organ = organ_label,
    Pollutant = exposure,
    `N candidates` = n_display,
    `Adverse events` = NA_integer_,
    `HR (95% CI)` = paste0("Mean difference ", result),
    `P value` = p_value_display,
    `Adjustment set` = adjustment_set
  )

tv_severity <- read_csv(
  file.path("output", "timevarying_pollution_severity_models_with_listing_year", "timevarying_pollution_severity_results.csv"),
  show_col_types = FALSE
) %>%
  mutate(
    analysis = "Time-updated organ severity association",
    model = "Pollution exposure as predictor of time-updated organ severity",
    result = format_est_ci(estimate, conf_low, conf_high, digits = 2),
    p_value_display = format_p(p_value)
  ) %>%
  transmute(
    Section = "Pollution-organ severity associations",
    Analysis = analysis,
    Model = model,
    Organ = organ_label,
    Pollutant = exposure,
    `N candidates` = candidate_organ_episodes,
    `Adverse events` = NA_integer_,
    `HR (95% CI)` = paste0("Mean difference ", result),
    `P value` = p_value_display,
    `Adjustment set` = adjustment_set
  )

supplement_table <- bind_rows(cox_table, baseline_severity, tv_severity)
write_csv(supplement_table, file.path(table_dir, "supplemental_updated_methods_sensitivity_table.csv"))
write_tsv(supplement_table, file.path(table_dir, "supplemental_updated_methods_sensitivity_table.tsv"))

notes <- c(
  "Supplemental table. Updated sensitivity analyses for prior-year baseline and time-varying Cox models.",
  "Baseline Cox models used 1-year exposure before listing and adjusted for age, sex, race, ZCTA-level SVI proxy, and transplant center strata unless listing center was omitted by design.",
  "Time-varying Cox models used monthly waitlist-period exposure intervals and adjusted for age, sex, race, ZCTA-level SVI proxy, organ-specific time-updated severity, and transplant center strata unless listing center was omitted by design.",
  "Multi-organ candidate sensitivity analyses either excluded candidates listed for more than one organ group or added an indicator for multi-organ candidate status.",
  "Time-varying multipollutant models were restricted to intervals with complete exposure data for included pollutants; daily PM2.5 and O3 availability limited the common multipollutant time-varying window to December 31, 2024.",
  "Pollution-organ severity association rows report adjusted mean differences in the organ-specific severity measure rather than hazard ratios."
)
writeLines(notes, file.path(table_dir, "supplemental_updated_methods_sensitivity_table_notes.txt"))

log_msg("Wrote updated methods sensitivity outputs to ", normalizePath(out_dir, winslash = "/"))
log_msg("Wrote supplemental sensitivity table to ", normalizePath(table_dir, winslash = "/"))
