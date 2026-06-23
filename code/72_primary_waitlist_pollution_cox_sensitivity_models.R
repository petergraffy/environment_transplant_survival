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
  library(readr)
  library(stringr)
  library(survival)
  library(tibble)
  library(tidyr)
})

analysis_path <- file.path(
  "output",
  "primary_waitlist_period_pollution_cox",
  "primary_waitlist_period_pollution_analysis_dataset.csv.gz"
)
pollution_dir <- file.path("data", "release", "air_pollution_zcta_parquet")
community_path <- file.path("data", "processed", "community", "zcta_acs_community_covariates_2005_2023.csv.gz")
out_dir <- file.path("output", "primary_waitlist_period_pollution_cox_sensitivity")
model_result_dir <- file.path(out_dir, "model_results")
dir.create(model_result_dir, recursive = TRUE, showWarnings = FALSE)

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
    select(zip, analysis_year, community_year, zcta_svi_proxy)
}

read_annual_pollution <- function() {
  pm25 <- read_parquet(file.path(pollution_dir, "air_pollution_zcta_pm25_monthly_2005_2023.parquet")) %>%
    transmute(zip = clean_zip(zip), year = as.integer(year), pm25_ug_m3) %>%
    group_by(zip, year) %>%
    summarise(pm25_n_months = sum(!is.na(pm25_ug_m3)), pm25_ug_m3 = mean(pm25_ug_m3, na.rm = TRUE), .groups = "drop") %>%
    filter(pm25_n_months == 12L, is.finite(pm25_ug_m3)) %>%
    select(zip, year, pm25_ug_m3)

  o3 <- read_parquet(file.path(pollution_dir, "air_pollution_zcta_o3_monthly_2005_2023.parquet")) %>%
    transmute(zip = clean_zip(zip), year = as.integer(year), o3_ppb) %>%
    group_by(zip, year) %>%
    summarise(o3_n_months = sum(!is.na(o3_ppb)), o3_ppb = mean(o3_ppb, na.rm = TRUE), .groups = "drop") %>%
    filter(o3_n_months == 12L, is.finite(o3_ppb)) %>%
    select(zip, year, o3_ppb)

  no2 <- read_parquet(file.path(pollution_dir, "air_pollution_zcta_no2_annual_2005_2025.parquet")) %>%
    transmute(zip = clean_zip(zip), year = as.integer(year), no2)

  list(pm25 = pm25, o3 = o3, no2 = no2)
}

primary_terms <- function(org, include_score = TRUE, include_center = TRUE, include_svi = FALSE,
                          include_kidney_clinical = TRUE) {
  terms <- c("age", "sex", "race", "listing_year")
  vars <- c("age", "sex", "race", "listing_year")

  if (org == "KI" && include_kidney_clinical) {
    terms <- c(terms, "kidney_no_dialysis_time", "kidney_dialysis_years", "kidney_diabetes")
    vars <- c(vars, "kidney_no_dialysis_time", "kidney_dialysis_years", "kidney_diabetes")
  } else if (include_score) {
    terms <- c(terms, "organ_score")
    vars <- c(vars, "organ_score")
  }

  if (include_svi) {
    terms <- c(terms, "zcta_svi_proxy")
    vars <- c(vars, "zcta_svi_proxy")
  }

  if (include_center) {
    terms <- c(terms, "strata(listing_center)")
    vars <- c(vars, "listing_center")
  }

  list(terms = terms, vars = vars)
}

format_adjustment <- function(org, include_score = TRUE, include_center = TRUE, include_svi = FALSE,
                              include_kidney_clinical = TRUE) {
  parts <- c("age", "sex", "race", "listing_year")
  if (org == "KI" && include_kidney_clinical) {
    parts <- c(parts, "no_dialysis_time", "dialysis_duration_years", "diabetes")
  } else if (include_score) {
    parts <- c(parts, "baseline_organ_score")
  }
  if (include_svi) parts <- c(parts, "zcta_svi_proxy")
  if (include_center) parts <- c(parts, "strata(listing_center)")
  paste(parts, collapse = " + ")
}

fit_sensitivity_cox <- function(dat, org, exposure_name, spec, sensitivity_name, sensitivity_label,
                                include_score = TRUE, include_center = TRUE, include_svi = FALSE,
                                include_kidney_clinical = TRUE) {
  adj <- primary_terms(
    org,
    include_score = include_score,
    include_center = include_center,
    include_svi = include_svi,
    include_kidney_clinical = include_kidney_clinical
  )
  vars_needed <- c(spec$followup, spec$event, adj$vars, spec$term)

  model_dat <- dat %>%
    filter(
      index_year <= spec$end_year,
      complete.cases(across(all_of(vars_needed))),
      .data[[spec$followup]] > 0
    ) %>%
    mutate(.followup_days = pmax(.data[[spec$followup]], 0.5), .event = .data[[spec$event]]) %>%
    droplevels()

  if (nrow(model_dat) == 0L || sum(model_dat$.event) < 50L) return(tibble())
  if (include_center && n_distinct(model_dat$listing_center) < 2L) return(tibble())

  form <- as.formula(paste("Surv(.followup_days, .event) ~", paste(c(spec$term, adj$terms), collapse = " + ")))
  log_msg(sensitivity_name, " ", org, " | ", exposure_name, " n=", nrow(model_dat), " adverse=", sum(model_dat$.event))
  fit <- coxph(form, data = model_dat, ties = "efron", x = FALSE, y = FALSE)

  tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(term == spec$term) %>%
    transmute(
      sensitivity = sensitivity_name,
      sensitivity_label = sensitivity_label,
      organ = org,
      organ_label = recode(org, !!!organ_labels),
      exposure = exposure_name,
      exposure_label = spec$label,
      exposure_window = spec$window,
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
      adjustment_set = format_adjustment(org, include_score, include_center, include_svi, include_kidney_clinical),
      variance_estimator = "model_based",
      center_adjustment = if_else(include_center, "stratified_baseline_hazard_by_listing_center", "none")
    )
}

fit_multipollutant_cox <- function(dat, org, exposure_names, specs, sensitivity_name, sensitivity_label,
                                   include_score = TRUE, include_center = TRUE, include_svi = FALSE,
                                   include_kidney_clinical = TRUE) {
  adj <- primary_terms(
    org,
    include_score = include_score,
    include_center = include_center,
    include_svi = include_svi,
    include_kidney_clinical = include_kidney_clinical
  )
  exposure_terms <- vapply(specs, `[[`, character(1), "term")
  names(exposure_terms) <- exposure_names
  exposure_lookup <- tibble(
    exposure = exposure_names,
    term = unname(exposure_terms),
    exposure_label = vapply(specs, `[[`, character(1), "label")
  )

  vars_needed <- c("followup_days_pm25", "event_pm25", adj$vars, unname(exposure_terms))

  model_dat <- dat %>%
    filter(
      index_year <= 2023L,
      complete.cases(across(all_of(vars_needed))),
      followup_days_pm25 > 0
    ) %>%
    mutate(.followup_days = pmax(followup_days_pm25, 0.5), .event = event_pm25) %>%
    droplevels()

  if (nrow(model_dat) == 0L || sum(model_dat$.event) < 50L) return(tibble())
  if (include_center && n_distinct(model_dat$listing_center) < 2L) return(tibble())

  form <- as.formula(paste(
    "Surv(.followup_days, .event) ~",
    paste(c(unname(exposure_terms), adj$terms), collapse = " + ")
  ))
  log_msg(
    sensitivity_name, " ", org, " | ",
    paste(exposure_names, collapse = "+"), " n=", nrow(model_dat),
    " adverse=", sum(model_dat$.event)
  )
  fit <- coxph(form, data = model_dat, ties = "efron", x = FALSE, y = FALSE)

  tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
    inner_join(exposure_lookup, by = "term") %>%
    transmute(
      sensitivity = sensitivity_name,
      sensitivity_label = sensitivity_label,
      organ = org,
      organ_label = recode(org, !!!organ_labels),
      exposure = exposure,
      exposure_label = exposure_label,
      exposure_window = "observed_waitlist_period_mean_multipollutant",
      exposure_data_end_year = 2023L,
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
      adjustment_set = format_adjustment(org, include_score, include_center, include_svi, include_kidney_clinical),
      variance_estimator = "model_based",
      center_adjustment = if_else(include_center, "stratified_baseline_hazard_by_listing_center", "none")
    )
}

log_msg("Reading primary analysis dataset")
analysis_dat <- read_csv(analysis_path, show_col_types = FALSE) %>%
  mutate(
    index_date = as.Date(index_date),
    observed_end_date = as.Date(observed_end_date),
    index_year = as.integer(format(index_date, "%Y")),
    listing_year = factor(as.integer(as.character(listing_year))),
    candidate_zip = clean_zip(candidate_zip),
    general_followup_days = as.numeric(observed_end_date - index_date),
    general_event = adverse_event,
    pm25_waitlist_5ug = pm25_waitlist_ug_m3 / 5,
    o3_waitlist_10ppb = o3_waitlist_ppb / 10,
    no2_waitlist_10unit = no2_waitlist / 10
  )

log_msg("Attaching ACS-derived ZCTA SVI proxy")
svi <- make_complete_acs_svi_proxy(community_path)
analysis_dat <- analysis_dat %>%
  left_join(
    svi,
    by = c("candidate_zip" = "zip", "index_year" = "analysis_year")
  )

log_msg("Attaching prior-year pollutant exposures")
annual_pollution <- read_annual_pollution()
prior_pm25 <- annual_pollution$pm25 %>%
  transmute(candidate_zip = zip, prior_exposure_year = year, pm25_prior_1y_5ug = pm25_ug_m3 / 5)
prior_o3 <- annual_pollution$o3 %>%
  transmute(candidate_zip = zip, prior_exposure_year = year, o3_prior_1y_10ppb = o3_ppb / 10)
prior_no2 <- annual_pollution$no2 %>%
  transmute(candidate_zip = zip, prior_exposure_year = year, no2_prior_1y_10unit = no2 / 10)

analysis_dat <- analysis_dat %>%
  mutate(prior_exposure_year = index_year - 1L) %>%
  left_join(prior_pm25, by = c("candidate_zip", "prior_exposure_year")) %>%
  left_join(prior_o3, by = c("candidate_zip", "prior_exposure_year")) %>%
  left_join(prior_no2, by = c("candidate_zip", "prior_exposure_year"))

waitlist_exposure_specs <- list(
  pm25 = list(term = "pm25_waitlist_5ug", label = "PM2.5 waitlist-period mean per 5 ug/m3", followup = "followup_days_pm25", event = "event_pm25", end_year = 2023L, window = "observed_waitlist_period_mean"),
  o3 = list(term = "o3_waitlist_10ppb", label = "Ozone waitlist-period mean per 10 ppb", followup = "followup_days_o3", event = "event_o3", end_year = 2023L, window = "observed_waitlist_period_mean"),
  no2 = list(term = "no2_waitlist_10unit", label = "NO2 waitlist-period mean per 10-unit increment", followup = "followup_days_no2", event = "event_no2", end_year = 2025L, window = "observed_waitlist_period_mean")
)

prior_exposure_specs <- list(
  pm25 = list(term = "pm25_prior_1y_5ug", label = "PM2.5 one-year prior mean per 5 ug/m3", followup = "general_followup_days", event = "general_event", end_year = 2024L, window = "one_year_prior_listing"),
  o3 = list(term = "o3_prior_1y_10ppb", label = "Ozone one-year prior mean per 10 ppb", followup = "general_followup_days", event = "general_event", end_year = 2024L, window = "one_year_prior_listing"),
  no2 = list(term = "no2_prior_1y_10unit", label = "NO2 one-year prior mean per 10-unit increment", followup = "general_followup_days", event = "general_event", end_year = 2025L, window = "one_year_prior_listing")
)

sensitivity_definitions <- tribble(
  ~sensitivity, ~sensitivity_label, ~spec_set, ~include_score, ~include_center, ~include_svi, ~include_kidney_clinical,
  "no_organ_score", "Primary waitlist-period exposure model without organ score", "waitlist", FALSE, TRUE, FALSE, TRUE,
  "no_listing_center", "Primary waitlist-period exposure model without listing-center strata", "waitlist", TRUE, FALSE, FALSE, TRUE,
  "primary_plus_svi", "Primary waitlist-period exposure model additionally adjusted for ACS-derived ZCTA SVI proxy", "waitlist", TRUE, TRUE, TRUE, TRUE,
  "prior_1y_no_score_no_center", "One-year prior exposure model without organ score, kidney clinical covariates, or listing-center strata", "prior", FALSE, FALSE, FALSE, FALSE
)

multipollutant_definitions <- tribble(
  ~sensitivity, ~sensitivity_label, ~exposure_names, ~include_score, ~include_center, ~include_svi, ~include_kidney_clinical,
  "multipollutant_pm25_no2", "Multipollutant waitlist-period model with PM2.5 and NO2", list(c("pm25", "no2")), TRUE, TRUE, FALSE, TRUE,
  "multipollutant_pm25_no2_o3", "Multipollutant waitlist-period model with PM2.5, NO2, and O3", list(c("pm25", "no2", "o3")), TRUE, TRUE, FALSE, TRUE
)

log_msg("Fitting sensitivity Cox models")
all_results <- list()
for (i in seq_len(nrow(sensitivity_definitions))) {
  def <- sensitivity_definitions[i, ]
  specs <- if (def$spec_set == "prior") prior_exposure_specs else waitlist_exposure_specs
  for (org in names(organ_labels)) {
    org_dat <- analysis_dat %>% filter(WL_ORG == org)
    for (exposure_name in names(specs)) {
      result <- fit_sensitivity_cox(
        org_dat,
        org,
        exposure_name,
        specs[[exposure_name]],
        def$sensitivity,
        def$sensitivity_label,
        include_score = def$include_score,
        include_center = def$include_center,
        include_svi = def$include_svi,
        include_kidney_clinical = def$include_kidney_clinical
      )
      result_path <- file.path(model_result_dir, paste0(def$sensitivity, "_", tolower(org), "_", exposure_name, ".csv"))
      write_csv(result, result_path)
      all_results[[length(all_results) + 1L]] <- result
      rm(result)
      invisible(gc())
    }
  }
}

for (i in seq_len(nrow(multipollutant_definitions))) {
  def <- multipollutant_definitions[i, ]
  for (org in names(organ_labels)) {
    org_dat <- analysis_dat %>% filter(WL_ORG == org)
    exposure_names <- unlist(def$exposure_names)
    specs <- waitlist_exposure_specs[exposure_names]
    result <- fit_multipollutant_cox(
      org_dat,
      org,
      exposure_names,
      specs,
      def$sensitivity,
      def$sensitivity_label,
      include_score = def$include_score,
      include_center = def$include_center,
      include_svi = def$include_svi,
      include_kidney_clinical = def$include_kidney_clinical
    )
    result_path <- file.path(
      model_result_dir,
      paste0(def$sensitivity, "_", tolower(org), ".csv")
    )
    write_csv(result, result_path)
    all_results[[length(all_results) + 1L]] <- result
    rm(result)
    invisible(gc())
  }
}

sensitivity_results <- bind_rows(all_results)
write_csv(sensitivity_results, file.path(out_dir, "primary_waitlist_period_pollution_cox_sensitivity_results.csv"))

write_csv(
  sensitivity_results %>%
    mutate(
      hazard_ratio_ci = sprintf("%.2f (%.2f to %.2f)", hazard_ratio, conf_low, conf_high),
      p_value_display = if_else(p_value < 0.001, "<.001", sprintf("%.3f", p_value))
    ) %>%
    select(sensitivity, organ_label, exposure, exposure_window, n, adverse_events, hazard_ratio_ci, p_value_display, adjustment_set),
  file.path(out_dir, "primary_waitlist_period_pollution_cox_sensitivity_table.csv")
)

log_msg("Wrote sensitivity Cox outputs to ", normalizePath(out_dir, winslash = "/"))
