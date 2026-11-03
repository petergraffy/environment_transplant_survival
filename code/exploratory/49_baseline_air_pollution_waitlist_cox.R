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
out_dir <- file.path("output", "baseline_air_pollution_waitlist_cox")
model_result_dir <- file.path(out_dir, "model_results")
dir.create(model_result_dir, recursive = TRUE, showWarnings = FALSE)

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

read_baseline_pollution <- function() {
  log_msg("Reading annual baseline air-pollution exposure tables")

  pm25 <- read_parquet(file.path(pollution_dir, "air_pollution_zcta_pm25_monthly_2005_2023.parquet")) %>%
    transmute(zip = clean_zip(zip), exposure_year = as.integer(year), pm25_ug_m3) %>%
    group_by(zip, exposure_year) %>%
    summarise(pm25_n_months = sum(!is.na(pm25_ug_m3)), pm25_ug_m3 = mean(pm25_ug_m3, na.rm = TRUE), .groups = "drop") %>%
    filter(pm25_n_months == 12L, is.finite(pm25_ug_m3))

  o3 <- read_parquet(file.path(pollution_dir, "air_pollution_zcta_o3_monthly_2005_2023.parquet")) %>%
    transmute(zip = clean_zip(zip), exposure_year = as.integer(year), o3_ppb) %>%
    group_by(zip, exposure_year) %>%
    summarise(o3_n_months = sum(!is.na(o3_ppb)), o3_ppb = mean(o3_ppb, na.rm = TRUE), .groups = "drop") %>%
    filter(o3_n_months == 12L, is.finite(o3_ppb))

  no2 <- read_parquet(file.path(pollution_dir, "air_pollution_zcta_no2_annual_2005_2025.parquet")) %>%
    transmute(zip = clean_zip(zip), exposure_year = as.integer(year), no2) %>%
    filter(!is.na(no2))

  list(pm25 = pm25, o3 = o3, no2 = no2)
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
    raw_event_date = coalesce(CAN_REM_DT, CAN_DEATH_DT, CAN_ENDWLFU),
    observed_end_date = if_else(is.na(raw_event_date) | raw_event_date > analysis_end_date, analysis_end_date, raw_event_date),
    event_observed = !is.na(raw_event_date) & raw_event_date <= analysis_end_date,
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
    followup_days = as.numeric(observed_end_date - index_date),
    age = as.numeric(CAN_AGE_AT_LISTING),
    sex = as_factor_missing(if_else(CAN_GENDER %in% c("M", "F"), CAN_GENDER, NA_character_)),
    race = as_factor_missing(CAN_RACE_SRTR),
    listing_center = as_factor_missing(CAN_LISTING_CTR_CD),
    listing_year = factor(index_year),
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
  mutate(
    followup_days = pmax(followup_days, 0.5),
    no2_exposure_year = index_year,
    pm25_exposure_year = index_year,
    o3_exposure_year = index_year,
    no2_prior_exposure_year = index_year - 1L,
    pm25_prior_exposure_year = index_year - 1L,
    o3_prior_exposure_year = index_year - 1L
  )

pollution <- read_baseline_pollution()

analysis_dat <- cohort %>%
  left_join(pollution$no2, by = c("candidate_zip" = "zip", "no2_exposure_year" = "exposure_year")) %>%
  rename(no2_listing_year = no2) %>%
  left_join(pollution$no2, by = c("candidate_zip" = "zip", "no2_prior_exposure_year" = "exposure_year")) %>%
  rename(no2_prior_year = no2) %>%
  left_join(pollution$pm25, by = c("candidate_zip" = "zip", "pm25_exposure_year" = "exposure_year")) %>%
  rename(pm25_listing_year_ug_m3 = pm25_ug_m3, pm25_listing_year_n_months = pm25_n_months) %>%
  left_join(pollution$pm25, by = c("candidate_zip" = "zip", "pm25_prior_exposure_year" = "exposure_year")) %>%
  rename(pm25_prior_year_ug_m3 = pm25_ug_m3, pm25_prior_year_n_months = pm25_n_months) %>%
  left_join(pollution$o3, by = c("candidate_zip" = "zip", "o3_exposure_year" = "exposure_year")) %>%
  rename(o3_listing_year_ppb = o3_ppb, o3_listing_year_n_months = o3_n_months) %>%
  left_join(pollution$o3, by = c("candidate_zip" = "zip", "o3_prior_exposure_year" = "exposure_year")) %>%
  rename(o3_prior_year_ppb = o3_ppb, o3_prior_year_n_months = o3_n_months) %>%
  mutate(
    no2_listing_prior_mean = rowMeans(cbind(no2_listing_year, no2_prior_year), na.rm = FALSE),
    pm25_listing_prior_mean_ug_m3 = rowMeans(cbind(pm25_listing_year_ug_m3, pm25_prior_year_ug_m3), na.rm = FALSE),
    o3_listing_prior_mean_ppb = rowMeans(cbind(o3_listing_year_ppb, o3_prior_year_ppb), na.rm = FALSE),
    no2_10unit = no2_listing_year / 10,
    pm25_5ug = pm25_listing_year_ug_m3 / 5,
    o3_10ppb = o3_listing_year_ppb / 10,
    no2_listing_prior_mean_10unit = no2_listing_prior_mean / 10,
    pm25_listing_prior_mean_5ug = pm25_listing_prior_mean_ug_m3 / 5,
    o3_listing_prior_mean_10ppb = o3_listing_prior_mean_ppb / 10
  )

write_csv(
  analysis_dat %>%
    select(
      waitlist_row_id, PERS_ID, PX_ID, WL_ORG, candidate_zip, index_date, observed_end_date,
      followup_days, event_type, adverse_event, transplant_or_improvement, other_exit,
      age, sex, race, listing_year, organ_score, organ_score_name, listing_center,
      no2_exposure_year, no2_listing_year, no2_prior_exposure_year, no2_prior_year, no2_listing_prior_mean,
      pm25_exposure_year, pm25_listing_year_ug_m3, pm25_listing_year_n_months,
      pm25_prior_exposure_year, pm25_prior_year_ug_m3, pm25_prior_year_n_months,
      pm25_listing_prior_mean_ug_m3,
      o3_exposure_year, o3_listing_year_ppb, o3_listing_year_n_months,
      o3_prior_exposure_year, o3_prior_year_ppb, o3_prior_year_n_months,
      o3_listing_prior_mean_ppb
    ),
  file.path(out_dir, "baseline_air_pollution_waitlist_analysis_dataset.csv.gz")
)

exposure_specs <- list(
  listing_year = list(
    no2 = list(term = "no2_10unit", label = "NO2 per 10-unit increment"),
    pm25 = list(term = "pm25_5ug", label = "PM2.5 per 5 ug/m3"),
    o3 = list(term = "o3_10ppb", label = "Ozone per 10 ppb")
  ),
  listing_prior_mean = list(
    no2 = list(term = "no2_listing_prior_mean_10unit", label = "NO2 listing/prior-year mean per 10-unit increment"),
    pm25 = list(term = "pm25_listing_prior_mean_5ug", label = "PM2.5 listing/prior-year mean per 5 ug/m3"),
    o3 = list(term = "o3_listing_prior_mean_10ppb", label = "Ozone listing/prior-year mean per 10 ppb")
  )
)

fit_baseline_cox <- function(dat, org, exposure_window, exposure_name, exposure_term, exposure_label, score_adjustment = TRUE) {
  vars_needed <- c("followup_days", "adverse_event", "age", "sex", "race", "listing_year", "listing_center", exposure_term)
  if (isTRUE(score_adjustment)) vars_needed <- c(vars_needed, "organ_score")
  model_dat <- dat %>%
    filter(complete.cases(across(all_of(vars_needed))), followup_days > 0) %>%
    droplevels()

  if (nrow(model_dat) == 0L || sum(model_dat$adverse_event) < 50L || n_distinct(model_dat$listing_center) < 2L) {
    return(tibble())
  }

  model_variant <- if (isTRUE(score_adjustment)) "with_organ_score" else "without_organ_score"
  log_msg("Baseline pollution Cox ", org, " | ", model_variant, " | ", exposure_window, " | ", exposure_name, " n=", nrow(model_dat), " adverse=", sum(model_dat$adverse_event))
  adjustment_terms <- c(exposure_term, "age", "sex", "race", "listing_year")
  if (isTRUE(score_adjustment)) adjustment_terms <- c(adjustment_terms, "organ_score")
  adjustment_terms <- c(adjustment_terms, "strata(listing_center)")
  form <- as.formula(paste(
    "Surv(followup_days, adverse_event) ~",
    paste(adjustment_terms, collapse = " + ")
  ))
  fit <- coxph(form, data = model_dat, ties = "efron", x = FALSE, y = FALSE)

  tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(term == exposure_term) %>%
    transmute(
      organ = org,
      organ_label = recode(org, !!!organ_labels),
      model_variant = model_variant,
      exposure_window = exposure_window,
      exposure = exposure_name,
      exposure_label = exposure_label,
      endpoint = "death_or_deterioration_delist_cause_specific",
      n = nrow(model_dat),
      people = n_distinct(model_dat$PERS_ID),
      centers = n_distinct(model_dat$listing_center),
      adverse_events = sum(model_dat$adverse_event),
      transplant_or_improvement_censors = sum(model_dat$event_type == "transplant_or_improvement"),
      other_exit_censors = sum(model_dat$event_type == "other_exit"),
      administrative_censors = sum(model_dat$event_type == "administrative_censor"),
      hazard_ratio = estimate,
      conf_low = conf.low,
      conf_high = conf.high,
      p_value = p.value,
      concordance = unname(summary(fit)$concordance[1]),
      adjustment_set = if_else(
        isTRUE(score_adjustment),
        "age + sex + race + listing_year + baseline_organ_score + strata(listing_center)",
        "age + sex + race + listing_year + strata(listing_center)"
      ),
      variance_estimator = "model_based",
      center_adjustment = "stratified_baseline_hazard_by_listing_center"
    )
}

log_msg("Fitting organ-specific baseline annual pollution cause-specific Cox models")
result_paths <- character()
for (org in target_organs) {
  org_dat <- analysis_dat %>% filter(WL_ORG == org)
  for (score_adjustment in c(TRUE, FALSE)) {
    model_variant <- if (isTRUE(score_adjustment)) "with_organ_score" else "without_organ_score"
    for (exposure_window in names(exposure_specs)) {
      for (exposure_name in names(exposure_specs[[exposure_window]])) {
        spec <- exposure_specs[[exposure_window]][[exposure_name]]
        result_path <- file.path(model_result_dir, paste0(tolower(org), "_", model_variant, "_", exposure_window, "_", exposure_name, ".csv"))
        result_paths <- c(result_paths, result_path)
        result <- fit_baseline_cox(org_dat, org, exposure_window, exposure_name, spec$term, spec$label, score_adjustment = score_adjustment)
        write_csv(result, result_path)
        rm(result)
        invisible(gc())
      }
    }
  }
}

baseline_results <- bind_rows(lapply(result_paths[file.exists(result_paths)], read_csv, show_col_types = FALSE))
write_csv(baseline_results, file.path(out_dir, "baseline_air_pollution_cox_results.csv"))

write_csv(
  analysis_dat %>%
    group_by(WL_ORG, organ_score_name) %>%
    summarise(
      n = n(),
      people = n_distinct(PERS_ID),
      centers = n_distinct(listing_center),
      adverse_events = sum(adverse_event),
      transplant_or_improvement_censors = sum(event_type == "transplant_or_improvement"),
      other_exit_censors = sum(event_type == "other_exit"),
      administrative_censors = sum(event_type == "administrative_censor"),
      median_followup_days = median(followup_days),
      complete_no2 = sum(!is.na(no2_10unit)),
      complete_pm25 = sum(!is.na(pm25_5ug)),
      complete_o3 = sum(!is.na(o3_10ppb)),
      complete_no2_listing_prior_mean = sum(!is.na(no2_listing_prior_mean_10unit)),
      complete_pm25_listing_prior_mean = sum(!is.na(pm25_listing_prior_mean_5ug)),
      complete_o3_listing_prior_mean = sum(!is.na(o3_listing_prior_mean_10ppb)),
      mean_no2_listing_year = mean(no2_listing_year, na.rm = TRUE),
      mean_pm25_listing_year = mean(pm25_listing_year_ug_m3, na.rm = TRUE),
      mean_o3_listing_year = mean(o3_listing_year_ppb, na.rm = TRUE),
      mean_no2_listing_prior = mean(no2_listing_prior_mean, na.rm = TRUE),
      mean_pm25_listing_prior = mean(pm25_listing_prior_mean_ug_m3, na.rm = TRUE),
      mean_o3_listing_prior = mean(o3_listing_prior_mean_ppb, na.rm = TRUE),
      .groups = "drop"
    ),
  file.path(out_dir, "baseline_air_pollution_cohort_summary.csv")
)

log_msg("Wrote baseline air-pollution waitlist Cox outputs to ", out_dir)
