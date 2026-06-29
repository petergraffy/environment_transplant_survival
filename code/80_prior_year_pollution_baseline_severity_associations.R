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
out_dir <- file.path("output", "prior_year_pollution_baseline_severity_associations")
cache_dir <- file.path(out_dir, "cache")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

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

parquet_files <- function(path) {
  list.files(path, pattern = "[.]parquet$", full.names = TRUE)
}

make_daily_annual <- function(path, value_col, out_col, cache_file) {
  if (file.exists(cache_file)) {
    return(read_csv(cache_file, show_col_types = FALSE) %>% mutate(zip = clean_zip(zip)))
  }
  open_dataset(parquet_files(path)) %>%
    transmute(zip = zip, year = year, value = .data[[value_col]]) %>%
    group_by(zip, year) %>%
    summarise(
      n_days = sum(!is.na(value)),
      value = mean(value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(n_days >= 300, is.finite(value)) %>%
    collect() %>%
    transmute(zip = clean_zip(zip), year = as.integer(year), !!out_col := value) %>%
    {
      write_csv(., cache_file)
      .
    }
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

read_prior_pollution <- function() {
  log_msg("Summarising daily PM2.5 and O3 to annual prior-year values")
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

attach_prior_pollution <- function(dat, prior_pollution) {
  dat %>%
    mutate(prior_exposure_year = listing_year_int - 1L) %>%
    left_join(prior_pollution$pm25, by = c("candidate_zip" = "zip", "prior_exposure_year" = "year")) %>%
    left_join(prior_pollution$o3, by = c("candidate_zip" = "zip", "prior_exposure_year" = "year")) %>%
    left_join(prior_pollution$no2, by = c("candidate_zip" = "zip", "prior_exposure_year" = "year")) %>%
    mutate(
      pm25_prior_5ug = pm25_prior_ug_m3 / 5,
      o3_prior_10ppb = o3_prior_ppb / 10,
      no2_prior_10ppb = no2_prior_ppb / 10,
      liver_meld_display = if_else(!is.na(liver_meld) & liver_meld > 1000, liver_meld - 6200, liver_meld)
    )
}

fit_severity_model <- function(dat, org, outcome, outcome_label, pollutant, exposure_term, exposure_label, family) {
  org_dat <- dat %>% filter(WL_ORG == org)
  rhs <- c(exposure_term, "age", "sex", "race", "listing_year", "zcta_svi_proxy", "listing_center")
  vars_needed <- c(outcome, rhs)
  model_dat <- org_dat %>%
    filter(complete.cases(across(all_of(vars_needed))))
  if (nrow(model_dat) == 0L) return(tibble())

  form <- as.formula(paste(outcome, "~", paste(rhs, collapse = " + ")))
  fit <- if (family == "binomial") {
    glm(form, data = model_dat, family = binomial())
  } else {
    lm(form, data = model_dat)
  }
  tidy(fit, conf.int = FALSE) %>%
    filter(term == exposure_term) %>%
    mutate(
      conf.low = estimate - 1.96 * std.error,
      conf.high = estimate + 1.96 * std.error,
      estimate_display = if_else(family == "binomial", exp(estimate), estimate),
      conf_low_display = if_else(family == "binomial", exp(conf.low), conf.low),
      conf_high_display = if_else(family == "binomial", exp(conf.high), conf.high)
    ) %>%
    transmute(
      organ = org,
      organ_label = recode(org, !!!organ_labels),
      outcome = outcome,
      outcome_label = outcome_label,
      outcome_family = family,
      pollutant = pollutant,
      exposure = exposure_label,
      n = nrow(model_dat),
      estimate = estimate_display,
      conf_low = conf_low_display,
      conf_high = conf_high_display,
      p_value = p.value,
      measure = if_else(family == "binomial", "odds_ratio", "mean_difference"),
      adjustment_set = "age + sex + race + listing_year + zcta_svi_proxy + listing_center"
    )
}

log_msg("Reading primary deduplicated cohort")
analysis_dat <- read_csv(analysis_path, show_col_types = FALSE) %>%
  mutate(
    candidate_zip = clean_zip(candidate_zip),
    index_date = as.Date(index_date),
    listing_year_int = as.integer(as.character(listing_year)),
    listing_year = factor(listing_year_int),
    sex = factor(sex),
    race = factor(race),
    listing_center = factor(listing_center)
  )

log_msg("Attaching ACS-derived ZCTA SVI proxy")
svi <- make_complete_acs_svi_proxy(community_path)
analysis_dat <- analysis_dat %>%
  left_join(svi, by = c("candidate_zip" = "zip", "listing_year_int" = "analysis_year"), suffix = c("", "_community")) %>%
  mutate(zcta_svi_proxy = coalesce(zcta_svi_proxy, zcta_svi_proxy_community)) %>%
  select(-any_of("zcta_svi_proxy_community"))

prior_pollution <- read_prior_pollution()
analysis_dat <- attach_prior_pollution(analysis_dat, prior_pollution)

severity_specs <- tribble(
  ~org, ~outcome, ~outcome_label, ~family,
  "KI", "kidney_dialysis_years", "Dialysis duration at listing, years", "gaussian",
  "KI", "kidney_diabetes", "Diabetes at listing", "binomial",
  "LI", "liver_meld_display", "MELD/PELD at listing", "gaussian",
  "HR", "organ_score", "US-candidate risk score proxy at listing", "gaussian",
  "LU", "organ_score", "Continuous lung waitlist urgency score proxy at listing", "gaussian"
)

pollutant_specs <- tribble(
  ~pollutant, ~term, ~label,
  "pm25", "pm25_prior_5ug", "Prior-year PM2.5 per 5 ug/m3",
  "o3", "o3_prior_10ppb", "Prior-year O3 per 10 ppb",
  "no2", "no2_prior_10ppb", "Prior-year NO2 per 10 ppb"
)

results <- bind_rows(lapply(seq_len(nrow(severity_specs)), function(i) {
  sev <- severity_specs[i, ]
  bind_rows(lapply(seq_len(nrow(pollutant_specs)), function(j) {
    pol <- pollutant_specs[j, ]
    fit_severity_model(
      analysis_dat,
      sev$org,
      sev$outcome,
      sev$outcome_label,
      pol$pollutant,
      pol$term,
      pol$label,
      sev$family
    )
  }))
}))

write_csv(results, file.path(out_dir, "prior_year_pollution_baseline_severity_associations.csv"))
write_csv(
  results %>%
    mutate(
      estimate_ci = if_else(
        measure == "odds_ratio",
        sprintf("OR %.2f (%.2f-%.2f)", estimate, conf_low, conf_high),
        sprintf("%.2f (%.2f-%.2f)", estimate, conf_low, conf_high)
      ),
      p_value_display = if_else(p_value < 0.001, "<.001", sprintf("%.3f", p_value))
    ) %>%
    select(organ_label, outcome_label, pollutant, exposure, n, measure, estimate_ci, p_value_display, adjustment_set),
  file.path(out_dir, "prior_year_pollution_baseline_severity_associations_table.csv")
)

log_msg("Wrote prior-year pollution baseline severity association outputs to ", normalizePath(out_dir, winslash = "/"))
