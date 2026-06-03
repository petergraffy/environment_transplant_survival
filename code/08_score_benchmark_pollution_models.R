#!/usr/bin/env Rscript

runtime_source <- file.path("code", "r_runtime.R")
if (!file.exists(runtime_source)) runtime_source <- "r_runtime.R"
if (file.exists(runtime_source)) {
  source(runtime_source)
  ensure_user_library()
}

suppressPackageStartupMessages({
  library(broom)
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

pubsaf_dir <- saf_paths$pubsaf_dir
supp_dir <- saf_paths$supp_dir
pollution_dir <- file.path("output", "zip_pollution")
community_dir <- file.path("data", "processed", "community")
out_dir <- file.path("output", "score_benchmark_pollution")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

analysis_start_year <- 2006L
analysis_end_year <- 2023L
target_organs <- c("HR", "KI", "LI", "LU")
horizon_days <- 42

log_msg <- function(...) {
  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
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

safe_log <- function(x) {
  log(pmax(as.numeric(x), 0.01))
}

safe_log1p <- function(x) {
  log(pmax(as.numeric(x), 0) + 1)
}

scale_01_to_pct <- function(x) {
  as.numeric(x) * 100
}

mode_value <- function(x) {
  ux <- na.omit(x)
  if (length(ux) == 0L) return(NA_character_)
  names(which.max(table(ux)))
}

read_pollution_annual <- function() {
  col_types <- cols(zip = col_character(), .default = col_guess())

  pm25 <- read_csv(file.path(pollution_dir, "all_pm25_zip.csv.gz"), col_types = col_types, show_col_types = FALSE) %>%
    group_by(zip = clean_zip(zip), year) %>%
    summarise(pm25_n_months = sum(!is.na(pm25_ug_m3)), pm25_ug_m3 = mean(pm25_ug_m3, na.rm = TRUE), .groups = "drop")

  o3 <- read_csv(file.path(pollution_dir, "all_o3_zip.csv.gz"), col_types = col_types, show_col_types = FALSE) %>%
    group_by(zip = clean_zip(zip), year) %>%
    summarise(o3_n_months = sum(!is.na(o3_ppb)), o3_ppb = mean(o3_ppb, na.rm = TRUE), .groups = "drop")

  no2 <- read_csv(file.path(pollution_dir, "all_no2_zip.csv.gz"), col_types = col_types, show_col_types = FALSE) %>%
    transmute(zip = clean_zip(zip), year, no2 = no2)

  pm25 %>%
    inner_join(o3, by = c("zip", "year")) %>%
    inner_join(no2, by = c("zip", "year")) %>%
    filter(pm25_n_months >= 12L, o3_n_months >= 12L, !is.na(pm25_ug_m3), !is.na(o3_ppb), !is.na(no2)) %>%
    mutate(
      pm25_5ug = pm25_ug_m3 / 5,
      o3_10ppb = o3_ppb / 10,
      no2_10unit = no2 / 10
    )
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

collapse_rare <- function(x, min_n = 50L) {
  y <- as.character(x)
  y[is.na(y) | y == ""] <- "Missing"
  tab <- table(y)
  y[!y %in% names(tab)[tab >= min_n]] <- "Other"
  factor(y)
}

median_impute <- function(x) {
  x <- as.numeric(x)
  med <- median(x, na.rm = TRUE)
  if (!is.finite(med)) med <- NA_real_
  if_else(is.na(x), med, x)
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
  "CAN_GENDER", "CAN_RACE_SRTR", "CAN_ETHNICITY_SRTR", "CAN_BMI", "CAN_ABO",
  "CAN_LISTING_CTR_CD", "CAN_INIT_STAT", "CAN_LAST_STAT", "CAN_MED_COND", "CAN_DGN",
  "CAN_ON_DIAL", "CAN_DIAL", "CAN_DIAL_DT", "CAN_GFR", "CAN_MOST_RECENT_CREAT",
  "CAN_INIT_SRTR_LAB_MELD", "CAN_LAST_SRTR_LAB_MELD", "CAN_TOT_BILI", "CAN_TOT_ALBUMIN",
  "CAN_LAST_SERUM_SODIUM", "CAN_FVC", "CAN_FEV1", "CAN_PCO2", "CAN_AT_REST_O2",
  "CAN_SIX_MIN_WALK", "CAN_PULM_ART_MEAN", "CAN_CARDIAC_OUTPUT", "CAN_VENTILATOR",
  "CAN_ON_VENTILATOR", "CAN_ECMO", "CAN_VAD_TAH", "CAN_VAD1", "CAN_VAD2",
  "CAN_CORTICOST_DEPND", "CAN_FUNCTN_STAT"
)

log_msg("Reading candidate ZIP supplement")
candidate_zip <- read_sas(file.path(supp_dir, "canzip2506.sas7bdat")) %>%
  transmute(PERS_ID, PX_ID, candidate_zip = clean_zip(CAN_PERM_ZIP)) %>%
  distinct(PERS_ID, PX_ID, .keep_all = TRUE)

log_msg("Reading pollution and community tables")
pollution_annual <- read_pollution_annual() %>%
  filter(year >= analysis_start_year - 1L, year <= analysis_end_year)

community <- read_csv(
  file.path(community_dir, "zcta_acs_community_covariates_2005_2023.csv.gz"),
  col_types = cols(zip = col_character(), .default = col_guess()),
  show_col_types = FALSE
) %>%
  transmute(
    zip = clean_zip(zip),
    analysis_year,
    median_household_income,
    pct_poverty = scale_01_to_pct(pct_poverty),
    pct_bachelor_plus = scale_01_to_pct(pct_bachelor_plus),
    pct_unemployed = scale_01_to_pct(pct_unemployed),
    pct_no_vehicle = scale_01_to_pct(pct_no_vehicle),
    pct_black = scale_01_to_pct(pct_black),
    pct_hispanic = scale_01_to_pct(pct_hispanic)
  )

log_msg("Reading candidate SAF files")
candidate_data <- candidate_files %>%
  rowwise() %>%
  mutate(data = list(
    read_sas(path, col_select = any_of(candidate_cols)) %>%
      mutate(candidate_group = candidate_group)
  )) %>%
  ungroup() %>%
  select(data) %>%
  tidyr::unnest(data)

log_msg("Reading heart status-justification US-CRS components")
heart_just <- read_sas(
  file.path(pubsaf_dir, "statjust_hr1a.sas7bdat"),
  col_select = any_of(c(
    "PX_ID", "WL_ORG", "CANHX_CHG_DT", "CANHX_STAT_CD", "CANHX_VAD", "CANHX_LVAD_TYPE",
    "CANHX_RVAD_TYPE", "CANHX_TAH_DT", "CANHX_ECMO", "CANHX_LAB_SERUM_CREAT",
    "CANHX_LAB_BILI", "CANHX_LAB_ALBUMIN", "CANHX_LAB_SODIUM", "CANHX_LAB_BNP"
  ))
) %>%
  filter(WL_ORG == "HR") %>%
  arrange(PX_ID, CANHX_CHG_DT) %>%
  group_by(PX_ID) %>%
  slice_tail(n = 1) %>%
  ungroup() %>%
  transmute(
    PX_ID,
    hr_stat_cd = CANHX_STAT_CD,
    hr_short_mcs = as.integer(flag_yes(CANHX_ECMO) == 1L | flag_yes(CANHX_RVAD_TYPE) == 1L),
    hr_durable_lvad = flag_yes(CANHX_LVAD_TYPE),
    hr_creatinine = CANHX_LAB_SERUM_CREAT,
    hr_bilirubin = CANHX_LAB_BILI,
    hr_albumin = CANHX_LAB_ALBUMIN,
    hr_sodium = CANHX_LAB_SODIUM,
    hr_bnp = CANHX_LAB_BNP
  )

baseline_exposure <- pollution_annual %>%
  transmute(
    candidate_zip = zip,
    prior_exposure_year = year,
    pm25_prior_1y_5ug = pm25_5ug,
    o3_prior_1y_10ppb = o3_10ppb,
    no2_prior_1y_10unit = no2_10unit
  )

current_exposure <- pollution_annual %>%
  transmute(
    candidate_zip = zip,
    index_year = year,
    pm25_current_5ug = pm25_5ug,
    o3_current_10ppb = o3_10ppb,
    no2_current_10unit = no2_10unit
  )

log_msg("Constructing benchmark cohort")
cohort <- candidate_data %>%
  filter(WL_ORG %in% target_organs) %>%
  left_join(candidate_zip, by = c("PERS_ID", "PX_ID")) %>%
  left_join(heart_just, by = "PX_ID") %>%
  mutate(
    index_date = coalesce(CAN_ACTIVATE_DT, CAN_LISTING_DT),
    index_year = as.integer(format(index_date, "%Y")),
    prior_exposure_year = index_year - 1L,
    event_date = coalesce(CAN_REM_DT, CAN_DEATH_DT, CAN_ENDWLFU),
    followup_days = as.numeric(event_date - index_date),
    adverse_waitlist_event = as.integer(CAN_REM_CD %in% c(8, 13)),
    death_event = as.integer(coalesce(CAN_REM_CD == 8, FALSE)),
    deterioration_event = as.integer(coalesce(CAN_REM_CD == 13, FALSE)),
    transplant_event = as.integer(CAN_REM_CD %in% c(4, 14, 15, 18, 19, 21, 22, 23)),
    age = CAN_AGE_AT_LISTING,
    sex = as_factor_missing(if_else(CAN_GENDER %in% c("M", "F"), CAN_GENDER, NA_character_)),
    race_srtr = as_factor_missing(CAN_RACE_SRTR),
    ethnicity_srtr = as_factor_missing(CAN_ETHNICITY_SRTR),
    listing_center = collapse_rare(CAN_LISTING_CTR_CD, min_n = 100L),
    index_year_factor = factor(index_year),
    init_status = collapse_rare(coalesce(as.character(hr_stat_cd), as.character(CAN_INIT_STAT)), min_n = 50L),
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
    score_name = case_when(
      WL_ORG == "HR" ~ "US-CRS proxy",
      WL_ORG == "KI" ~ "Dialysis + age",
      WL_ORG == "LI" ~ "MELD",
      WL_ORG == "LU" ~ "LAS/CAS component proxy"
    ),
    benchmark_score = case_when(
      WL_ORG == "HR" ~ us_crs_proxy,
      WL_ORG == "KI" ~ kidney_dialysis_age_score,
      WL_ORG == "LI" ~ liver_meld,
      WL_ORG == "LU" ~ lung_las_cas_component_score
    )
  ) %>%
  filter(
    index_year >= analysis_start_year,
    index_year <= analysis_end_year,
    CAN_SOURCE %in% c("A", "R"),
    !is.na(index_date),
    !is.na(event_date),
    !is.na(followup_days),
    followup_days >= 0,
    !is.na(candidate_zip)
  ) %>%
  mutate(followup_days = pmax(followup_days, 0.5)) %>%
  left_join(baseline_exposure, by = c("candidate_zip", "prior_exposure_year")) %>%
  left_join(current_exposure, by = c("candidate_zip", "index_year")) %>%
  left_join(community, by = c("candidate_zip" = "zip", "index_year" = "analysis_year"))

community_terms <- c(
  "ns(median_household_income, df = 3)", "ns(pct_poverty, df = 3)",
  "ns(pct_bachelor_plus, df = 3)", "ns(pct_unemployed, df = 3)",
  "ns(pct_no_vehicle, df = 3)", "ns(pct_black, df = 3)", "ns(pct_hispanic, df = 3)"
)

pollution_prior_terms <- c("pm25_prior_1y_5ug", "o3_prior_1y_10ppb", "no2_prior_1y_10unit")
pollution_current_terms <- c("pm25_current_5ug", "o3_current_10ppb", "no2_current_10unit")

write_csv(
  cohort %>%
    group_by(WL_ORG, score_name) %>%
    summarise(
      n = n(),
      people = n_distinct(PERS_ID),
      adverse_events = sum(adverse_waitlist_event),
      deaths = sum(death_event),
      deteriorations = sum(deterioration_event),
      complete_score = sum(!is.na(benchmark_score)),
      complete_prior_pollution = sum(complete.cases(across(all_of(pollution_prior_terms)))),
      complete_current_pollution = sum(complete.cases(across(all_of(pollution_current_terms)))),
      .groups = "drop"
    ),
  file.path(out_dir, "score_benchmark_cohort_summary.csv")
)

extract_formula_vars <- function(terms) {
  vars <- unique(unlist(str_extract_all(terms, "[A-Za-z][A-Za-z0-9_]*")))
  setdiff(vars, c("ns", "df", "strata", "cluster"))
}

fit_model <- function(org, model_name, terms) {
  dat <- cohort %>% filter(WL_ORG == org)
  vars_needed <- intersect(extract_formula_vars(terms), names(dat))
  dat <- dat %>% filter(complete.cases(across(all_of(c(vars_needed, "benchmark_score"))))) %>% droplevels()
  if (nrow(dat) == 0L || sum(dat$adverse_waitlist_event) < 100L) return(NULL)

  rhs <- paste(terms, collapse = " + ")
  form <- as.formula(paste("Surv(followup_days, adverse_waitlist_event) ~", rhs))
  log_msg("Fitting ", org, " ", model_name, " n=", nrow(dat), " adverse=", sum(dat$adverse_waitlist_event))
  fit <- coxph(form, data = dat, ties = "efron", robust = TRUE, x = FALSE, y = FALSE)
  s <- summary(fit)

  tibble(
    organ = org,
    score_name = unique(dat$score_name),
    model = model_name,
    n = nrow(dat),
    adverse_events = sum(dat$adverse_waitlist_event),
    concordance = unname(s$concordance[1]),
    loglik = unname(fit$loglik[2]),
    aic = AIC(fit)
  ) %>%
    bind_cols(
      tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
        filter(term %in% c("benchmark_score", pollution_prior_terms, pollution_current_terms)) %>%
        select(term, estimate, conf.low, conf.high, p.value) %>%
        tidyr::nest(coef = everything())
    )
}

fit_horizon_model <- function(org, model_name, terms, horizon = horizon_days) {
  dat <- cohort %>%
    filter(WL_ORG == org) %>%
    mutate(
      horizon_followup_days = pmin(followup_days, horizon),
      horizon_adverse_event = as.integer(adverse_waitlist_event == 1L & followup_days <= horizon)
    )

  vars_needed <- intersect(extract_formula_vars(terms), names(dat))
  dat <- dat %>%
    filter(complete.cases(across(all_of(c(vars_needed, "benchmark_score"))))) %>%
    droplevels()
  if (nrow(dat) == 0L || sum(dat$horizon_adverse_event) < 50L) return(NULL)

  rhs <- paste(terms, collapse = " + ")
  form <- as.formula(paste("Surv(horizon_followup_days, horizon_adverse_event) ~", rhs))
  log_msg("Fitting ", org, " ", horizon, "-day ", model_name, " n=", nrow(dat), " adverse=", sum(dat$horizon_adverse_event))
  fit <- coxph(form, data = dat, ties = "efron", robust = TRUE, x = FALSE, y = FALSE)
  s <- summary(fit)

  tibble(
    organ = org,
    score_name = unique(dat$score_name),
    endpoint = paste0(horizon, "_day_adverse_waitlist"),
    model = model_name,
    n = nrow(dat),
    adverse_events = sum(dat$horizon_adverse_event),
    concordance = unname(s$concordance[1]),
    loglik = unname(fit$loglik[2]),
    aic = AIC(fit)
  ) %>%
    bind_cols(
      tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
        filter(term %in% c("benchmark_score", pollution_prior_terms, pollution_current_terms)) %>%
        select(term, estimate, conf.low, conf.high, p.value) %>%
        tidyr::nest(coef = everything())
    )
}

base_terms <- c("benchmark_score")
year_terms <- c("benchmark_score", "strata(index_year_factor)", "cluster(PERS_ID)")
year_community_terms <- c("benchmark_score", "strata(index_year_factor)", community_terms, "cluster(PERS_ID)")
year_center_terms <- c("benchmark_score", "strata(listing_center)", "strata(index_year_factor)", "cluster(PERS_ID)")
context_terms <- c("benchmark_score", "strata(listing_center)", "strata(index_year_factor)", community_terms, "cluster(PERS_ID)")
score_prior_terms <- c(base_terms, pollution_prior_terms, "cluster(PERS_ID)")
year_prior_terms <- c(year_terms, pollution_prior_terms)
year_community_prior_terms <- c(year_community_terms, pollution_prior_terms)
year_center_prior_terms <- c(year_center_terms, pollution_prior_terms)
prior_terms <- c(context_terms, pollution_prior_terms)
prior_current_terms <- c(context_terms, pollution_prior_terms, pollution_current_terms)

model_grid <- list(
  score_only = base_terms,
  score_prior_pollution = score_prior_terms,
  score_year_prior_pollution = year_prior_terms,
  score_year_community_prior_pollution = year_community_prior_terms,
  score_year_center_prior_pollution = year_center_prior_terms,
  score_center_year_community = context_terms,
  score_center_year_community_prior_pollution = prior_terms,
  score_center_year_community_prior_current_pollution = prior_current_terms
)

results <- bind_rows(lapply(target_organs, function(org) {
  bind_rows(lapply(names(model_grid), function(model_name) fit_model(org, model_name, model_grid[[model_name]])))
}))

model_summary <- results %>%
  select(-coef)

coef_results <- results %>%
  select(organ, score_name, model, n, adverse_events, concordance, coef) %>%
  tidyr::unnest(coef)

write_csv(model_summary, file.path(out_dir, "score_benchmark_model_performance.csv"))
write_csv(coef_results, file.path(out_dir, "score_benchmark_pollution_coefficients.csv"))

write_csv(
  model_summary %>%
    group_by(organ) %>%
    arrange(match(model, names(model_grid)), .by_group = TRUE) %>%
    mutate(
      delta_concordance_vs_score_only = concordance - first(concordance),
      delta_aic_vs_score_only = aic - first(aic)
    ) %>%
    ungroup(),
  file.path(out_dir, "score_benchmark_incremental_performance.csv")
)

heart_horizon_results <- bind_rows(lapply(names(model_grid), function(model_name) {
  fit_horizon_model("HR", model_name, model_grid[[model_name]], horizon = horizon_days)
}))

heart_horizon_summary <- heart_horizon_results %>% select(-coef)
heart_horizon_coef <- heart_horizon_results %>%
  select(organ, score_name, endpoint, model, n, adverse_events, concordance, coef) %>%
  tidyr::unnest(coef)

write_csv(heart_horizon_summary, file.path(out_dir, "heart_uscrs_42day_model_performance.csv"))
write_csv(heart_horizon_coef, file.path(out_dir, "heart_uscrs_42day_pollution_coefficients.csv"))
write_csv(
  heart_horizon_summary %>%
    arrange(match(model, names(model_grid))) %>%
    mutate(
      delta_concordance_vs_score_only = concordance - first(concordance),
      delta_aic_vs_score_only = aic - first(aic)
    ),
  file.path(out_dir, "heart_uscrs_42day_incremental_performance.csv")
)

message("Wrote benchmark score/pollution comparison outputs to ", out_dir)
