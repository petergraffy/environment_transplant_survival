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
community_path <- file.path("data", "processed", "community", "zcta_acs_community_covariates_2005_2023.csv.gz")
release_dir <- file.path("data", "release")
annual_pollution_dir <- file.path(release_dir, "air_pollution_zcta_parquet")
prior_cache_dir <- file.path("output", "prior_year_pollution_cox_svi", "cache")
out_dir <- file.path("output", "prior_year_initial_listing_cox_model_set")
model_dir <- file.path(out_dir, "model_results")
table_dir <- file.path("output", "tables", "prior_year_initial_listing_cox_model_set")
dir.create(prior_cache_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

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

format_p <- function(p) {
  case_when(
    is.na(p) ~ "",
    p < 0.001 ~ "<.001",
    TRUE ~ sprintf("%.3f", p)
  )
}

format_hr <- function(hr, low, high) {
  sprintf("%.2f (%.2f-%.2f)", hr, low, high)
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



severity_terms_for_organ <- function(org) {
  if (org == "KI") {
    return(c("kidney_no_dialysis_time", "kidney_dialysis_years", "kidney_diabetes"))
  }
  "organ_score"
}

model_terms_for_tier <- function(org, model_tier) {
  if (model_tier == "Unadjusted") return(character())

  terms <- c("age", "sex", "race", "zcta_svi_proxy", "strata(listing_center)")
  if (model_tier == "Adjusted + SDOH + organ score at listing") {
    terms <- c(terms, severity_terms_for_organ(org))
  }
  terms
}

model_vars_for_tier <- function(org, model_tier) {
  vars <- model_terms_for_tier(org, model_tier)
  vars <- vars[!grepl("^strata\\(", vars)]
  if ("strata(listing_center)" %in% model_terms_for_tier(org, model_tier)) {
    vars <- c(vars, "listing_center")
  }
  unique(vars)
}

fit_initial_listing_model <- function(dat, org, model_tier, model_type, exposure_terms, exposure_labels) {
  adjustment_terms <- model_terms_for_tier(org, model_tier)
  adjustment_vars <- model_vars_for_tier(org, model_tier)
  vars_needed <- c("followup_days", "adverse_event", exposure_terms, adjustment_vars)

  model_dat <- dat %>%
    filter(
      WL_ORG == org,
      complete.cases(across(all_of(vars_needed))),
      followup_days > 0
    ) %>%
    mutate(.followup_days = pmax(followup_days, 0.5)) %>%
    droplevels()

  if (nrow(model_dat) == 0L || sum(model_dat$adverse_event, na.rm = TRUE) < 50L) {
    return(tibble())
  }
  if ("listing_center" %in% vars_needed && n_distinct(model_dat$listing_center) < 2L) {
    return(tibble())
  }

  log_msg(
    "Initial-listing prior-year Cox ", org, " | ", model_tier, " | ",
    model_type, " n=", nrow(model_dat), " adverse=", sum(model_dat$adverse_event)
  )

  rhs_terms <- c(exposure_terms, adjustment_terms)
  form <- as.formula(paste(
    "Surv(.followup_days, adverse_event) ~",
    paste(rhs_terms, collapse = " + ")
  ))

  fit <- coxph(form, data = model_dat, ties = "efron", x = FALSE, y = FALSE)

  tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(term %in% exposure_terms) %>%
    transmute(
      analysis = "Initial-listing prior-year cause-specific Cox",
      model_tier = model_tier,
      model_type = model_type,
      organ = org,
      organ_label = recode(org, !!!organ_labels),
      exposure = term,
      exposure_label = unname(exposure_labels[term]),
      exposure_window = "previous_365_days_pm25_o3_or_12_complete_months_no2",
      endpoint = "death_or_deterioration_delist_cause_specific",
      n = nrow(model_dat),
      people = n_distinct(model_dat$PERS_ID),
      centers = if ("listing_center" %in% vars_needed) n_distinct(model_dat$listing_center) else NA_integer_,
      adverse_events = sum(model_dat$adverse_event),
      hazard_ratio = estimate,
      conf_low = conf.low,
      conf_high = conf.high,
      p_value = p.value,
      concordance = unname(summary(fit)$concordance[1]),
      adjustment_set = if (length(adjustment_terms)) paste(adjustment_terms, collapse = " + ") else "None",
      variance_estimator = "model_based",
      center_adjustment = if ("listing_center" %in% vars_needed) "stratified_baseline_hazard_by_listing_center" else "none"
    )
}

fit_initial_listing_model_coefficients <- function(dat, org, model_tier, model_type, exposure_terms, exposure_labels) {
  adjustment_terms <- model_terms_for_tier(org, model_tier)
  adjustment_vars <- model_vars_for_tier(org, model_tier)
  vars_needed <- c("followup_days", "adverse_event", exposure_terms, adjustment_vars)

  model_dat <- dat %>%
    filter(
      WL_ORG == org,
      complete.cases(across(all_of(vars_needed))),
      followup_days > 0
    ) %>%
    mutate(.followup_days = pmax(followup_days, 0.5)) %>%
    droplevels()

  if (nrow(model_dat) == 0L || sum(model_dat$adverse_event, na.rm = TRUE) < 50L) {
    return(tibble())
  }
  if ("listing_center" %in% vars_needed && n_distinct(model_dat$listing_center) < 2L) {
    return(tibble())
  }

  rhs_terms <- c(exposure_terms, adjustment_terms)
  form <- as.formula(paste(
    "Surv(.followup_days, adverse_event) ~",
    paste(rhs_terms, collapse = " + ")
  ))

  fit <- coxph(form, data = model_dat, ties = "efron", x = FALSE, y = FALSE)

  tidy(fit, exponentiate = FALSE, conf.int = FALSE) %>%
    mutate(
      conf_low = estimate - 1.96 * std.error,
      conf_high = estimate + 1.96 * std.error,
      hazard_ratio = exp(estimate),
      hr_conf_low = exp(conf_low),
      hr_conf_high = exp(conf_high),
      term_type = case_when(
        term %in% exposure_terms ~ "pollutant",
        term == "age" ~ "covariate",
        startsWith(term, "sex") ~ "covariate",
        startsWith(term, "race") ~ "covariate",
        term == "zcta_svi_proxy" ~ "covariate",
        term %in% severity_terms_for_organ(org) ~ "organ_severity",
        TRUE ~ "covariate"
      ),
      exposure_label = if_else(term %in% names(exposure_labels), unname(exposure_labels[term]), NA_character_)
    ) %>%
    transmute(
      analysis = "Initial-listing prior-year cause-specific Cox",
      model_tier = model_tier,
      model_type = model_type,
      organ = org,
      organ_label = recode(org, !!!organ_labels),
      term,
      term_type,
      exposure_label,
      exposure_window = "previous_365_days_pm25_o3_or_12_complete_months_no2",
      endpoint = "death_or_deterioration_delist_cause_specific",
      n = nrow(model_dat),
      people = n_distinct(model_dat$PERS_ID),
      centers = if ("listing_center" %in% vars_needed) n_distinct(model_dat$listing_center) else NA_integer_,
      adverse_events = sum(model_dat$adverse_event),
      coefficient = estimate,
      std_error = std.error,
      conf_low,
      conf_high,
      hazard_ratio,
      hr_conf_low,
      hr_conf_high,
      p_value = p.value,
      concordance = unname(summary(fit)$concordance[1]),
      adjustment_set = if (length(adjustment_terms)) paste(adjustment_terms, collapse = " + ") else "None",
      variance_estimator = "model_based",
      center_adjustment = if ("listing_center" %in% vars_needed) "stratified_baseline_hazard_by_listing_center" else "none"
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
    listing_center = factor(listing_center),
    kidney_no_dialysis_time = as.integer(kidney_no_dialysis_time),
    kidney_diabetes = as.integer(kidney_diabetes)
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

exposure_term_labels <- c(
  pm25_prior_5ug = "Prior-year PM2.5 per 5 ug/m3",
  o3_prior_10ppb = "Prior-year O3 per 10 ppb",
  no2_prior_10ppb = "Prior-year NO2 per 10 ppb"
)

model_tiers <- c(
  "Unadjusted",
  "Adjusted + SDOH",
  "Adjusted + SDOH + organ score at listing"
)
model_specs <- list(
  "Single-pollutant PM2.5" = "pm25_prior_5ug",
  "Single-pollutant O3" = "o3_prior_10ppb",
  "Single-pollutant NO2" = "no2_prior_10ppb",
  "Multipollutant PM2.5 + NO2" = c("pm25_prior_5ug", "no2_prior_10ppb"),
  "Multipollutant PM2.5 + NO2 + O3" = c("pm25_prior_5ug", "no2_prior_10ppb", "o3_prior_10ppb")
)

results <- list()
coefficient_results <- list()
for (org in target_organs) {
  for (tier in model_tiers) {
    for (model_type in names(model_specs)) {
      terms <- model_specs[[model_type]]
      result <- fit_initial_listing_model(
        analysis_dat,
        org,
        tier,
        model_type,
        terms,
        exposure_term_labels[terms]
      )
      result_path <- file.path(
        model_dir,
        paste0(
          tolower(org), "_",
          str_replace_all(str_to_lower(tier), "[^a-z0-9]+", "_"), "_",
          str_replace_all(str_to_lower(model_type), "[^a-z0-9]+", "_"),
          ".csv"
        )
      )
      write_csv(result, result_path)
      results[[length(results) + 1L]] <- result

      coefficient_results[[length(coefficient_results) + 1L]] <- fit_initial_listing_model_coefficients(
        analysis_dat,
        org,
        tier,
        model_type,
        terms,
        exposure_term_labels[terms]
      )
    }
  }
}

results <- bind_rows(results)
coefficient_results <- bind_rows(coefficient_results)
write_csv(results, file.path(out_dir, "prior_year_initial_listing_cox_model_set_results.csv"))
write_csv(coefficient_results, file.path(out_dir, "prior_year_initial_listing_cox_model_set_full_coefficients.csv"))

wide_table <- results %>%
  mutate(
    result = format_hr(hazard_ratio, conf_low, conf_high),
    p_value_display = format_p(p_value),
    pollutant = case_when(
      exposure == "pm25_prior_5ug" ~ "PM2.5",
      exposure == "o3_prior_10ppb" ~ "O3",
      exposure == "no2_prior_10ppb" ~ "NO2",
      TRUE ~ exposure
    ),
    model_column = paste(model_tier, model_type, sep = " | ")
  ) %>%
  select(organ_label, pollutant, model_column, result) %>%
  distinct() %>%
  pivot_wider(names_from = model_column, values_from = result, values_fill = "") %>%
  arrange(factor(organ_label, levels = unname(organ_labels)), factor(pollutant, levels = c("PM2.5", "O3", "NO2"))) %>%
  rename(Organ = organ_label, Pollutant = pollutant)

write_csv(wide_table, file.path(table_dir, "prior_year_initial_listing_cox_model_set_wide_table.csv"))
write_tsv(wide_table, file.path(table_dir, "prior_year_initial_listing_cox_model_set_wide_table.tsv"))

long_table <- results %>%
  mutate(
    `HR (95% CI)` = format_hr(hazard_ratio, conf_low, conf_high),
    `P value` = format_p(p_value)
  ) %>%
  transmute(
    Analysis = analysis,
    `Model tier` = model_tier,
    `Model type` = model_type,
    Organ = organ_label,
    Pollutant = exposure_label,
    `N candidates` = n,
    `Adverse events` = adverse_events,
    `HR (95% CI)` = `HR (95% CI)`,
    `P value` = `P value`,
    `Adjustment set` = adjustment_set
  )
write_csv(long_table, file.path(table_dir, "prior_year_initial_listing_cox_model_set_long_table.csv"))
write_tsv(long_table, file.path(table_dir, "prior_year_initial_listing_cox_model_set_long_table.tsv"))

coefficient_table <- coefficient_results %>%
  mutate(
    `Coefficient (95% CI)` = sprintf("%.4f (%.4f to %.4f)", coefficient, conf_low, conf_high),
    `HR (95% CI)` = format_hr(hazard_ratio, hr_conf_low, hr_conf_high),
    `P value` = format_p(p_value)
  ) %>%
  transmute(
    Analysis = analysis,
    `Model tier` = model_tier,
    `Model type` = model_type,
    Organ = organ_label,
    Term = term,
    `Term type` = term_type,
    `Exposure label` = coalesce(exposure_label, ""),
    `N candidates` = n,
    `Adverse events` = adverse_events,
    `Coefficient (95% CI)` = `Coefficient (95% CI)`,
    `HR (95% CI)` = `HR (95% CI)`,
    `P value` = `P value`,
    `Adjustment set` = adjustment_set
  )
write_csv(coefficient_table, file.path(table_dir, "prior_year_initial_listing_cox_model_set_full_coefficients.csv"))
write_tsv(coefficient_table, file.path(table_dir, "prior_year_initial_listing_cox_model_set_full_coefficients.tsv"))

cohort_summary <- bind_rows(lapply(names(model_specs), function(model_type) {
  terms <- model_specs[[model_type]]
  bind_rows(lapply(model_tiers, function(tier) {
    bind_rows(lapply(target_organs, function(org) {
      vars_needed <- c("followup_days", "adverse_event", terms, model_vars_for_tier(org, tier))
      analysis_dat %>%
        filter(
          WL_ORG == org,
          complete.cases(across(all_of(vars_needed))),
          followup_days > 0
        ) %>%
        summarise(
          model_tier = tier,
          model_type = model_type,
          organ = org,
          organ_label = recode(org, !!!organ_labels),
          n = n(),
          people = n_distinct(PERS_ID),
          adverse_events = sum(adverse_event, na.rm = TRUE),
          median_followup_days = median(followup_days, na.rm = TRUE),
          .groups = "drop"
        )
    }))
  }))
}))
write_csv(cohort_summary, file.path(out_dir, "prior_year_initial_listing_cox_model_set_cohort_summary.csv"))

notes <- c(
  "Initial-listing prior-year cause-specific Cox model set.",
  "All models use mean PM2.5/O3 during the 365 days before initial listing, or day-weighted NO2 during the 12 complete calendar months before the listing month. Annual NO2 approximates months before monthly coverage begins in 2019.",
  "Unadjusted models include pollutant exposure term(s) only.",
  "Adjusted + SDOH models include age, sex, race, ZCTA-level SVI proxy, and transplant center strata.",
  "Adjusted + SDOH + organ score at listing models additionally include baseline organ_score for heart, liver, and lung candidates; kidney models additionally include no dialysis time, dialysis duration in years, and diabetes.",
  "Single-pollutant models were fit separately for PM2.5, O3, and NO2. Multipollutant models included PM2.5 + NO2 and PM2.5 + NO2 + O3.",
  "Values are hazard ratios with 95% confidence intervals; PM2.5 is per 5 ug/m3, O3 per 10 ppb, and NO2 per 10 ppb."
)
writeLines(notes, file.path(table_dir, "prior_year_initial_listing_cox_model_set_notes.txt"))

log_msg("Wrote initial-listing prior-year Cox model set to ", normalizePath(out_dir, winslash = "/"))
log_msg("Wrote tables to ", normalizePath(table_dir, winslash = "/"))
