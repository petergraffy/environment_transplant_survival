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
out_dir <- file.path("output", "pared_primary_etiologic")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

analysis_start_year <- 2006L
analysis_end_year <- 2023L
target_organs <- c("HR", "KI", "LI", "LU")
pollutant_terms <- c(
  pm25_prior_1y_5ug = "PM2.5 per 5 ug/m3",
  o3_prior_1y_10ppb = "Ozone per 10 ppb",
  no2_prior_1y_10unit = "NO2 per 10-unit increase"
)

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

scale_01_to_pct <- function(x) {
  x <- as.numeric(x)
  if_else(!is.na(x) & x <= 1, x * 100, x)
}

collapse_rare <- function(x, min_n = 50L) {
  y <- as.character(x)
  y[is.na(y) | y == ""] <- "Missing"
  tab <- table(y)
  y[!y %in% names(tab)[tab >= min_n]] <- "Other"
  factor(y)
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
    mutate(pm25_5ug = pm25_ug_m3 / 5, o3_10ppb = o3_ppb / 10, no2_10unit = no2 / 10)
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
  "CAN_GENDER", "CAN_RACE_SRTR", "CAN_ETHNICITY_SRTR", "CAN_LISTING_CTR_CD",
  "CAN_INIT_STAT", "CAN_ON_DIAL", "CAN_DIAL", "CAN_DIAL_DT", "CAN_MOST_RECENT_CREAT",
  "CAN_INIT_SRTR_LAB_MELD", "CAN_LAST_SRTR_LAB_MELD", "CAN_TOT_BILI", "CAN_TOT_ALBUMIN",
  "CAN_LAST_SERUM_SODIUM", "CAN_FVC", "CAN_FEV1", "CAN_PCO2", "CAN_AT_REST_O2",
  "CAN_SIX_MIN_WALK", "CAN_PULM_ART_MEAN", "CAN_CARDIAC_OUTPUT", "CAN_VENTILATOR",
  "CAN_ON_VENTILATOR", "CAN_ECMO", "CAN_VAD_TAH", "CAN_VAD1", "CAN_VAD2",
  "CAN_CORTICOST_DEPND"
)

log_msg("Reading ZIP, pollution, and SAF inputs")
candidate_zip <- read_sas(file.path(supp_dir, "canzip2506.sas7bdat")) %>%
  transmute(PERS_ID, PX_ID, candidate_zip = clean_zip(CAN_PERM_ZIP)) %>%
  distinct(PERS_ID, PX_ID, .keep_all = TRUE)

pollution_prior <- read_pollution_annual() %>%
  filter(year >= analysis_start_year - 1L, year <= analysis_end_year - 1L) %>%
  transmute(
    candidate_zip = zip,
    prior_exposure_year = year,
    pm25_prior_1y_5ug = pm25_5ug,
    o3_prior_1y_10ppb = o3_10ppb,
    no2_prior_1y_10unit = no2_10unit
  )

center_region <- read_sas(file.path(pubsaf_dir, "institution.sas7bdat"), col_select = any_of(c("CTR_CD", "REGION"))) %>%
  transmute(
    listing_center_raw = as.character(CTR_CD),
    optn_region = factor(REGION)
  ) %>%
  filter(!is.na(listing_center_raw), !is.na(optn_region)) %>%
  distinct(listing_center_raw, .keep_all = TRUE)

community_covariates <- read_csv(
  file.path(community_dir, "zcta_acs_community_covariates_2005_2023.csv.gz"),
  col_types = cols(zip = col_character(), .default = col_guess()),
  show_col_types = FALSE
) %>%
  transmute(
    candidate_zip = clean_zip(zip),
    index_year = as.integer(analysis_year),
    median_household_income = as.numeric(median_household_income),
    pct_poverty = scale_01_to_pct(pct_poverty),
    pct_bachelor_plus = scale_01_to_pct(pct_bachelor_plus),
    pct_unemployed = scale_01_to_pct(pct_unemployed),
    pct_no_vehicle = scale_01_to_pct(pct_no_vehicle),
    pct_black = scale_01_to_pct(pct_black),
    pct_hispanic = scale_01_to_pct(pct_hispanic)
  ) %>%
  distinct(candidate_zip, index_year, .keep_all = TRUE)

candidate_data <- candidate_files %>%
  rowwise() %>%
  mutate(data = list(read_sas(path, col_select = any_of(candidate_cols)) %>% mutate(candidate_group = candidate_group))) %>%
  ungroup() %>%
  select(data) %>%
  tidyr::unnest(data)

heart_just <- read_sas(
  file.path(pubsaf_dir, "statjust_hr1a.sas7bdat"),
  col_select = any_of(c(
    "PX_ID", "WL_ORG", "CANHX_CHG_DT", "CANHX_STAT_CD", "CANHX_LVAD_TYPE",
    "CANHX_RVAD_TYPE", "CANHX_ECMO", "CANHX_LAB_SERUM_CREAT", "CANHX_LAB_BILI",
    "CANHX_LAB_ALBUMIN", "CANHX_LAB_SODIUM", "CANHX_LAB_BNP"
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

log_msg("Constructing pared primary cohort")
cohort <- candidate_data %>%
  filter(WL_ORG %in% target_organs) %>%
  left_join(candidate_zip, by = c("PERS_ID", "PX_ID")) %>%
  left_join(heart_just, by = "PX_ID") %>%
  mutate(
    index_date = coalesce(CAN_ACTIVATE_DT, CAN_LISTING_DT),
    index_year = as.integer(format(index_date, "%Y")),
    index_year_centered = index_year - 2015,
    prior_exposure_year = index_year - 1L,
    event_date = coalesce(CAN_REM_DT, CAN_DEATH_DT, CAN_ENDWLFU),
    followup_days = as.numeric(event_date - index_date),
    adverse_waitlist_event = as.integer(CAN_REM_CD %in% c(8, 13)),
    transplant_event = as.integer(CAN_REM_CD %in% c(4, 14, 15, 18, 19, 21, 22, 23)),
    age = CAN_AGE_AT_LISTING,
    sex = as_factor_missing(if_else(CAN_GENDER %in% c("M", "F"), CAN_GENDER, NA_character_)),
    race = collapse_rare(CAN_RACE_SRTR, min_n = 100L),
    listing_center_raw = as.character(CAN_LISTING_CTR_CD),
    listing_center = collapse_rare(CAN_LISTING_CTR_CD, min_n = 100L),
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
      WL_ORG == "HR" ~ "US-CRS proxy",
      WL_ORG == "KI" ~ "Dialysis + age",
      WL_ORG == "LI" ~ "MELD",
      WL_ORG == "LU" ~ "LAS/CAS component proxy"
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
  left_join(pollution_prior, by = c("candidate_zip", "prior_exposure_year")) %>%
  left_join(center_region, by = "listing_center_raw") %>%
  left_join(community_covariates, by = c("candidate_zip", "index_year")) %>%
  mutate(income_10k = median_household_income / 10000) %>%
  filter(complete.cases(
    age, sex, race, index_year_centered, organ_score, listing_center, optn_region,
    income_10k, pct_poverty, pct_bachelor_plus, pct_unemployed, pct_no_vehicle, pct_black, pct_hispanic,
    pm25_prior_1y_5ug, o3_prior_1y_10ppb, no2_prior_1y_10unit
  ))

write_csv(
  cohort %>%
    group_by(WL_ORG, organ_score_name) %>%
    summarise(
      n = n(),
      people = n_distinct(PERS_ID),
      centers = n_distinct(listing_center_raw),
      optn_regions = n_distinct(optn_region),
      adverse_events = sum(adverse_waitlist_event),
      transplants = sum(transplant_event),
      median_followup_days = median(followup_days),
      .groups = "drop"
    ),
  file.path(out_dir, "pared_primary_cohort_summary.csv")
)

community_terms <- c(
  "income_10k",
  "pct_poverty",
  "pct_bachelor_plus",
  "pct_unemployed",
  "pct_no_vehicle",
  "pct_black",
  "pct_hispanic"
)

stage_terms <- list(
  pollutant_only = character(),
  demographics_year = c("age", "sex", "race", "index_year_centered"),
  clinical_score = c("age", "sex", "race", "index_year_centered", "organ_score"),
  optn_region = c("age", "sex", "race", "index_year_centered", "organ_score", "optn_region"),
  community = c("age", "sex", "race", "index_year_centered", "organ_score", "optn_region", community_terms),
  center_frailty_sensitivity = c(
    "age", "sex", "race", "index_year_centered", "organ_score", "optn_region", community_terms,
    "frailty(listing_center, distribution = 'gaussian')"
  )
)

fit_stage <- function(dat, organ, pollutant_term, pollutant_label, stage_name, terms) {
  rhs <- c(pollutant_term, terms)
  form <- as.formula(paste("Surv(followup_days, adverse_waitlist_event) ~", paste(rhs, collapse = " + ")))
  log_msg("Fitting ", organ, " | ", pollutant_label, " | ", stage_name, " n=", nrow(dat), " adverse=", sum(dat$adverse_waitlist_event))
  has_frailty <- str_detect(stage_name, "frailty")
  fit <- coxph(form, data = dat, ties = "efron", robust = !has_frailty, x = FALSE, y = FALSE)
  s <- summary(fit)
  tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(term == pollutant_term) %>%
    transmute(
      organ,
      organ_score_name = unique(dat$organ_score_name),
      pollutant = pollutant_label,
      pollutant_term,
      stage = stage_name,
      n = nrow(dat),
      adverse_events = sum(dat$adverse_waitlist_event),
      hazard_ratio = estimate,
      conf_low = conf.low,
      conf_high = conf.high,
      p_value = p.value,
      concordance = unname(s$concordance[1]),
      aic = AIC(fit)
    )
}

cv_stage <- function(dat, organ, pollutant_term, pollutant_label, stage_name, terms, k = 5L) {
  set.seed(20260529 + match(organ, target_organs))
  people <- unique(dat$PERS_ID)
  fold_map <- tibble(PERS_ID = people, fold = sample(rep(seq_len(k), length.out = length(people))))
  dat <- dat %>% left_join(fold_map, by = "PERS_ID")
  rhs <- c(pollutant_term, terms)
  form <- as.formula(paste("Surv(followup_days, adverse_waitlist_event) ~", paste(rhs, collapse = " + ")))
  fold_results <- bind_rows(lapply(seq_len(k), function(fold_id) {
    train <- dat %>% filter(fold != fold_id) %>% droplevels()
    test <- dat %>% filter(fold == fold_id) %>% droplevels()
    fit <- coxph(form, data = train, ties = "efron", x = FALSE, y = FALSE)
    test$risk_score <- -as.numeric(predict(fit, newdata = test, type = "lp"))
    cidx <- concordance(Surv(followup_days, adverse_waitlist_event) ~ risk_score, data = test)$concordance
    tibble(
      organ,
      pollutant = pollutant_label,
      pollutant_term,
      stage = stage_name,
      fold = fold_id,
      n_train = nrow(train),
      n_test = nrow(test),
      adverse_test = sum(test$adverse_waitlist_event),
      cv_concordance = cidx
    )
  }))
  fold_results
}

fixed_stage_names <- setdiff(names(stage_terms), "center_frailty_sensitivity")

model_results <- bind_rows(lapply(target_organs, function(org) {
  dat <- cohort %>% filter(WL_ORG == org) %>% droplevels()
  bind_rows(lapply(names(pollutant_terms), function(pollutant_term) {
    bind_rows(lapply(names(stage_terms), function(stage_name) {
      fit_stage(dat, org, pollutant_term, pollutant_terms[[pollutant_term]], stage_name, stage_terms[[stage_name]])
    }))
  }))
}))

cv_results <- bind_rows(lapply(target_organs, function(org) {
  dat <- cohort %>% filter(WL_ORG == org) %>% droplevels()
  bind_rows(lapply(names(pollutant_terms), function(pollutant_term) {
    bind_rows(lapply(fixed_stage_names, function(stage_name) {
      cv_stage(dat, org, pollutant_term, pollutant_terms[[pollutant_term]], stage_name, stage_terms[[stage_name]])
    }))
  }))
}))

cv_summary <- cv_results %>%
  group_by(organ, pollutant, pollutant_term, stage) %>%
  summarise(
    mean_cv_concordance = mean(cv_concordance, na.rm = TRUE),
    sd_cv_concordance = sd(cv_concordance, na.rm = TRUE),
    total_test_adverse = sum(adverse_test),
    .groups = "drop"
  )

clean_staged_table <- model_results %>%
  left_join(cv_summary, by = c("organ", "pollutant", "pollutant_term", "stage")) %>%
  arrange(organ, pollutant, match(stage, names(stage_terms)))

write_csv(model_results, file.path(out_dir, "pared_primary_model_pollutant_coefficients.csv"))
write_csv(cv_results, file.path(out_dir, "pared_primary_model_cv_folds.csv"))
write_csv(cv_summary, file.path(out_dir, "pared_primary_model_cv_summary.csv"))
write_csv(clean_staged_table, file.path(out_dir, "pared_primary_clean_staged_table.csv"))

message("Wrote pared primary etiologic model outputs to ", out_dir)
