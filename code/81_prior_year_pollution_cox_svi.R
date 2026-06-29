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
  library(dplyr)
  library(readr)
  library(stringr)
  library(survival)
  library(tibble)
})

analysis_path <- file.path(
  "output",
  "primary_waitlist_period_pollution_cox",
  "primary_waitlist_period_pollution_analysis_dataset.csv.gz"
)
community_path <- file.path("data", "processed", "community", "zcta_acs_community_covariates_2005_2023.csv.gz")
release_dir <- file.path("data", "release")
annual_pollution_dir <- file.path(release_dir, "air_pollution_zcta_parquet")
out_dir <- file.path("output", "prior_year_pollution_cox_svi")
model_result_dir <- file.path(out_dir, "model_results")
cache_dir <- file.path(out_dir, "cache")
dir.create(model_result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

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

parquet_files <- function(path) {
  list.files(path, pattern = "[.]parquet$", full.names = TRUE)
}

make_complete_acs_svi_proxy <- function(path) {
  community <- read_csv(path, show_col_types = FALSE) %>%
    mutate(
      zip = clean_zip(zip),
      analysis_year = as.integer(analysis_year)
    )

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

make_daily_annual <- function(path, value_col, out_col, cache_file) {
  if (file.exists(cache_file)) {
    return(read_csv(cache_file, show_col_types = FALSE) %>% mutate(zip = clean_zip(zip)))
  }
  annual <- open_dataset(parquet_files(path)) %>%
    transmute(zip = zip, year = year, value = .data[[value_col]]) %>%
    group_by(zip, year) %>%
    summarise(
      n_days = sum(!is.na(value)),
      value = mean(value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(n_days >= 300, is.finite(value)) %>%
    collect() %>%
    transmute(zip = clean_zip(zip), year = as.integer(year), !!out_col := value)
  write_csv(annual, cache_file)
  annual
}

read_prior_pollution <- function() {
  pm25 <- make_daily_annual(
    file.path(release_dir, "lghap_pm25_zcta_daily_parquet"),
    "pm25_ug_m3",
    "pm25_prior_ug_m3",
    file.path(cache_dir, "pm25_daily_annual_zcta.csv.gz")
  )
  o3 <- make_daily_annual(
    file.path(release_dir, "o3_zcta_daily_parquet"),
    "o3_ppb",
    "o3_prior_ppb",
    file.path(cache_dir, "o3_daily_annual_zcta.csv.gz")
  )
  no2 <- read_parquet(file.path(annual_pollution_dir, "air_pollution_zcta_no2_annual_2005_2025.parquet")) %>%
    transmute(zip = clean_zip(zip), year = as.integer(year), no2_prior_ppb = no2)
  list(pm25 = pm25, o3 = o3, no2 = no2)
}

fit_prior_year_cox <- function(dat, org, pollutant, exposure_term, exposure_label) {
  vars_needed <- c("followup_days", "adverse_event", exposure_term, "age", "sex", "race", "zcta_svi_proxy", "listing_center")
  model_dat <- dat %>%
    filter(
      WL_ORG == org,
      complete.cases(across(all_of(vars_needed))),
      followup_days > 0
    ) %>%
    mutate(.followup_days = pmax(followup_days, 0.5)) %>%
    droplevels()

  if (nrow(model_dat) == 0L || sum(model_dat$adverse_event) < 50L || n_distinct(model_dat$listing_center) < 2L) {
    return(tibble())
  }

  log_msg("Prior-year pollution Cox ", org, " | ", pollutant, " n=", nrow(model_dat), " adverse=", sum(model_dat$adverse_event))
  form <- as.formula(paste(
    "Surv(.followup_days, adverse_event) ~",
    paste(c(exposure_term, "age", "sex", "race", "zcta_svi_proxy", "strata(listing_center)"), collapse = " + ")
  ))
  fit <- coxph(form, data = model_dat, ties = "efron", x = FALSE, y = FALSE)

  tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(term == exposure_term) %>%
    transmute(
      organ = org,
      organ_label = recode(org, !!!organ_labels),
      exposure = pollutant,
      exposure_label = exposure_label,
      exposure_window = "calendar_year_before_listing",
      endpoint = "death_or_deterioration_delist_cause_specific",
      n = nrow(model_dat),
      people = n_distinct(model_dat$PERS_ID),
      centers = n_distinct(model_dat$listing_center),
      adverse_events = sum(model_dat$adverse_event),
      hazard_ratio = estimate,
      conf_low = conf.low,
      conf_high = conf.high,
      p_value = p.value,
      concordance = unname(summary(fit)$concordance[1]),
      adjustment_set = "age + sex + race + zcta_svi_proxy + strata(listing_center)",
      variance_estimator = "model_based",
      center_adjustment = "stratified_baseline_hazard_by_listing_center"
    )
}

log_msg("Reading primary deduplicated cohort")
analysis_dat <- read_csv(analysis_path, show_col_types = FALSE) %>%
  mutate(
    candidate_zip = clean_zip(candidate_zip),
    index_date = as.Date(index_date),
    observed_end_date = as.Date(observed_end_date),
    listing_year_int = as.integer(as.character(listing_year)),
    prior_exposure_year = listing_year_int - 1L,
    followup_days = as.numeric(observed_end_date - index_date),
    sex = factor(sex),
    race = factor(race),
    listing_center = factor(listing_center)
  ) %>%
  filter(WL_ORG %in% target_organs)

log_msg("Attaching ACS-derived ZCTA SVI proxy")
svi <- make_complete_acs_svi_proxy(community_path)
analysis_dat <- analysis_dat %>%
  left_join(svi, by = c("candidate_zip" = "zip", "listing_year_int" = "analysis_year"), suffix = c("", "_community")) %>%
  mutate(zcta_svi_proxy = coalesce(zcta_svi_proxy, zcta_svi_proxy_community)) %>%
  select(-any_of("zcta_svi_proxy_community"))

log_msg("Attaching prior-year pollution")
prior_pollution <- read_prior_pollution()
analysis_dat <- analysis_dat %>%
  left_join(prior_pollution$pm25, by = c("candidate_zip" = "zip", "prior_exposure_year" = "year")) %>%
  left_join(prior_pollution$o3, by = c("candidate_zip" = "zip", "prior_exposure_year" = "year")) %>%
  left_join(prior_pollution$no2, by = c("candidate_zip" = "zip", "prior_exposure_year" = "year")) %>%
  mutate(
    pm25_prior_5ug = pm25_prior_ug_m3 / 5,
    o3_prior_10ppb = o3_prior_ppb / 10,
    no2_prior_10ppb = no2_prior_ppb / 10
  )

pollutant_specs <- tribble(
  ~pollutant, ~term, ~label,
  "pm25", "pm25_prior_5ug", "Prior-year PM2.5 per 5 ug/m3",
  "o3", "o3_prior_10ppb", "Prior-year O3 per 10 ppb",
  "no2", "no2_prior_10ppb", "Prior-year NO2 per 10 ppb"
)

result_paths <- character()
for (org in target_organs) {
  for (i in seq_len(nrow(pollutant_specs))) {
    spec <- pollutant_specs[i, ]
    result <- fit_prior_year_cox(analysis_dat, org, spec$pollutant, spec$term, spec$label)
    result_path <- file.path(model_result_dir, paste0(tolower(org), "_", spec$pollutant, "_prior_year_cox.csv"))
    write_csv(result, result_path)
    result_paths <- c(result_paths, result_path)
  }
}

results <- bind_rows(lapply(result_paths[file.exists(result_paths)], read_csv, show_col_types = FALSE))
write_csv(results, file.path(out_dir, "prior_year_pollution_cox_svi_results.csv"))
write_csv(
  results %>%
    mutate(
      hazard_ratio_ci = sprintf("%.2f (%.2f-%.2f)", hazard_ratio, conf_low, conf_high),
      p_value_display = if_else(p_value < 0.001, "<.001", sprintf("%.3f", p_value))
    ) %>%
    select(organ_label, exposure, exposure_label, n, adverse_events, hazard_ratio_ci, p_value_display, adjustment_set),
  file.path(out_dir, "prior_year_pollution_cox_svi_table.csv")
)

write_csv(
  bind_rows(lapply(seq_len(nrow(pollutant_specs)), function(i) {
    spec <- pollutant_specs[i, ]
    analysis_dat %>%
      filter(!is.na(.data[[spec$term]]), complete.cases(age, sex, race, zcta_svi_proxy, listing_center), followup_days > 0) %>%
      group_by(WL_ORG) %>%
      summarise(
        exposure = spec$pollutant,
        n = n(),
        people = n_distinct(PERS_ID),
        adverse_events = sum(adverse_event, na.rm = TRUE),
        median_followup_days = median(followup_days, na.rm = TRUE),
        mean_prior_exposure = mean(.data[[spec$term]], na.rm = TRUE),
        .groups = "drop"
      )
  })),
  file.path(out_dir, "prior_year_pollution_cohort_summary.csv")
)

log_msg("Wrote prior-year pollution Cox outputs to ", normalizePath(out_dir, winslash = "/"))
