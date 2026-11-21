#!/usr/bin/env Rscript

runtime_source <- file.path("code", "r_runtime.R")
if (!file.exists(runtime_source)) runtime_source <- "r_runtime.R"
if (file.exists(runtime_source)) {
  source(runtime_source)
  ensure_user_library()
}

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(haven)
  library(readr)
  library(stringr)
  library(tibble)
  library(tidyr)
})

saf_paths_source <- file.path("code", "saf_paths.R")
if (!file.exists(saf_paths_source)) saf_paths_source <- "saf_paths.R"
source(saf_paths_source)
saf_paths <- get_saf_paths(release = "q1_2026")
assert_saf_files(saf_paths, include_stathist = TRUE)

pubsaf_dir <- saf_paths$pubsaf_dir
out_dir <- file.path("output", "primary_waitlist_period_pollution_cox")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

analysis_start_year <- as.integer(Sys.getenv("ENV_ANALYSIS_START_YEAR", "2005"))
pollution_end_year <- c(pm25 = 2023L, o3 = 2023L, no2 = 2025L)
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

durable_vad_brand_codes <- c(
  202L, 205L, 206L, 207L, 208L, 209L, 210L, 212L, 213L, 214L,
  223L, 224L, 233L, 236L, 239L, 240L, 312L, 313L, 314L, 315L,
  316L, 319L, 322L, 327L, 330L, 333L, 334L
)

flag_durable_vad_brand <- function(...) {
  vals <- list(...)
  Reduce(`|`, lapply(vals, function(x) suppressWarnings(as.integer(as.character(x))) %in% durable_vad_brand_codes)) %>%
    as.integer()
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

assign_overlap_episode <- function(index_date, observed_end_date) {
  episode_id <- integer(length(index_date))
  if (!length(index_date)) return(episode_id)
  current_episode <- 1L
  current_end <- observed_end_date[[1]]
  episode_id[[1]] <- current_episode
  if (length(index_date) == 1L) return(episode_id)
  for (i in seq.int(2L, length(index_date))) {
    if (is.na(index_date[[i]]) || is.na(current_end) || index_date[[i]] > current_end) {
      current_episode <- current_episode + 1L
      current_end <- observed_end_date[[i]]
    } else {
      current_end <- max(current_end, observed_end_date[[i]], na.rm = TRUE)
    }
    episode_id[[i]] <- current_episode
  }
  episode_id
}

collapse_overlap_episode <- function(episode_dat) {
  episode_dat <- episode_dat %>% arrange(index_date, observed_end_date, waitlist_row_id)
  baseline <- episode_dat %>% slice_head(n = 1)

  adverse_dates <- episode_dat$observed_end_date[episode_dat$adverse_event == 1L]
  good_dates <- episode_dat$observed_end_date[episode_dat$transplant_or_improvement == 1L]
  terminal_dates <- c(adverse_dates, good_dates)

  if (length(terminal_dates) > 0L) {
    episode_end <- min(terminal_dates, na.rm = TRUE)
    adverse_event <- as.integer(any(episode_dat$adverse_event == 1L & episode_dat$observed_end_date == episode_end, na.rm = TRUE))
    transplant_or_improvement <- as.integer(adverse_event == 0L && any(episode_dat$transplant_or_improvement == 1L & episode_dat$observed_end_date == episode_end, na.rm = TRUE))
    other_exit <- 0L
  } else {
    episode_end <- max(episode_dat$observed_end_date, na.rm = TRUE)
    adverse_event <- 0L
    transplant_or_improvement <- 0L
    other_exit <- as.integer(any(episode_dat$other_exit == 1L, na.rm = TRUE))
  }

  baseline %>%
    mutate(
      observed_end_date = episode_end,
      raw_event_date = if_else(adverse_event == 1L | transplant_or_improvement == 1L | other_exit == 1L, episode_end, as.Date(NA)),
      event_observed = adverse_event == 1L | transplant_or_improvement == 1L | other_exit == 1L,
      adverse_event = adverse_event,
      transplant_or_improvement = transplant_or_improvement,
      other_exit = other_exit,
      event_type = factor(
        case_when(
          adverse_event == 1L ~ "death_or_deterioration_delist",
          transplant_or_improvement == 1L ~ "transplant_or_improvement",
          other_exit == 1L ~ "other_exit",
          TRUE ~ "administrative_censor"
        ),
        levels = c("administrative_censor", "death_or_deterioration_delist", "transplant_or_improvement", "other_exit")
      ),
      episode_listing_count = nrow(episode_dat),
      episode_center_count = n_distinct(episode_dat$listing_center),
      episode_px_ids = paste(unique(episode_dat$PX_ID), collapse = ";")
    )
}

event_priority <- function(adverse_event, transplant_or_improvement, other_exit) {
  case_when(
    adverse_event == 1L ~ 3L,
    transplant_or_improvement == 1L ~ 2L,
    other_exit == 1L ~ 1L,
    TRUE ~ 0L
  )
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
  "CAN_GENDER", "CAN_RACE_SRTR", "CAN_LISTING_CTR_CD", "CAN_ON_DIAL",
  "CAN_DIAL", "CAN_DIAL_DT", "CAN_DIAB", "CAN_DIAB_TY", "CAN_MOST_RECENT_CREAT", "CAN_INIT_SRTR_LAB_MELD",
  "CAN_LAST_SRTR_LAB_MELD", "CAN_TOT_BILI", "CAN_TOT_ALBUMIN",
  "CAN_LAST_SERUM_SODIUM", "CAN_FVC", "CAN_FEV1", "CAN_PCO2", "CAN_AT_REST_O2",
  "CAN_SIX_MIN_WALK", "CAN_PULM_ART_MEAN", "CAN_CARDIAC_OUTPUT",
  "CAN_VENTILATOR", "CAN_ON_VENTILATOR", "CAN_ECMO", "CAN_VAD1",
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
    "PX_ID", "WL_ORG", "CANHX_CHG_DT", "CANHX_RVAD_TYPE",
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
    hr_creatinine = CANHX_LAB_SERUM_CREAT,
    hr_bilirubin = CANHX_LAB_BILI,
    hr_albumin = CANHX_LAB_ALBUMIN,
    hr_sodium = CANHX_LAB_SODIUM,
    hr_bnp = CANHX_LAB_BNP
  )

heart_riskstrat <- read_sas(
  file.path(saf_paths$supp_dir, "Thoracic Registration", "RiskStratDataHR.sas7bdat"),
  col_select = any_of(c(
    "px_id", "RiskStratId", "ChangeDate", "HrSevFailCreatinine",
    "HrSevFailBilirubin", "HrSevFailAlbumin", "HrSevFailSodium",
    "HrSevFailBnp", "HrSevFailNtBnpType"
  ))
) %>%
  rename(PX_ID = px_id) %>%
  arrange(PX_ID, ChangeDate, RiskStratId) %>%
  group_by(PX_ID) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  transmute(
    PX_ID,
    hr_riskstrat_creatinine = HrSevFailCreatinine,
    hr_riskstrat_bilirubin = HrSevFailBilirubin,
    hr_riskstrat_albumin = HrSevFailAlbumin,
    hr_riskstrat_sodium = HrSevFailSodium,
    hr_riskstrat_bnp = HrSevFailBnp,
    hr_riskstrat_bnp_type = HrSevFailNtBnpType
  )

log_msg("Constructing waitlist cohort")
max_end_date <- as.Date(sprintf("%04d-12-31", max(pollution_end_year)))
cohort <- candidate_data %>%
  filter(WL_ORG %in% target_organs) %>%
  left_join(candidate_zip, by = c("PERS_ID", "PX_ID")) %>%
  left_join(heart_baseline, by = "PX_ID") %>%
  left_join(heart_riskstrat, by = "PX_ID") %>%
  mutate(
    waitlist_row_id = row_number(),
    index_date = CAN_LISTING_DT,
    index_year = as.integer(format(index_date, "%Y")),
    raw_event_date = coalesce(CAN_REM_DT, CAN_DEATH_DT, CAN_ENDWLFU),
    observed_end_date = if_else(is.na(raw_event_date) | raw_event_date > max_end_date, max_end_date, raw_event_date),
    event_observed = !is.na(raw_event_date) & raw_event_date <= max_end_date,
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
    age = as.numeric(CAN_AGE_AT_LISTING),
    sex = as_factor_missing(if_else(CAN_GENDER %in% c("M", "F"), CAN_GENDER, NA_character_)),
    race = as_factor_missing(CAN_RACE_SRTR),
    listing_center = as_factor_missing(CAN_LISTING_CTR_CD),
    listing_year = factor(index_year),
    kidney_on_dialysis = flag_yes(coalesce(as.character(CAN_ON_DIAL), as.character(CAN_DIAL))),
    kidney_diabetes = as.integer(
      CAN_DIAB %in% c(2, 3, 4) |
        CAN_DIAB_TY %in% c(2, 3, 4, 5)
    ),
    dialysis_years_raw = as.numeric(index_date - CAN_DIAL_DT) / 365.25,
    dialysis_years_raw = if_else(is.na(dialysis_years_raw) | dialysis_years_raw < 0, NA_real_, dialysis_years_raw),
    kidney_no_dialysis_time = as.integer(kidney_on_dialysis != 1L | is.na(dialysis_years_raw) | dialysis_years_raw == 0),
    kidney_dialysis_years = if_else(kidney_no_dialysis_time == 1L, 0, dialysis_years_raw),
    dialysis_years = kidney_dialysis_years,
    liver_meld = coalesce(CAN_INIT_SRTR_LAB_MELD, CAN_LAST_SRTR_LAB_MELD),
    hr_creatinine_for_score = coalesce(hr_riskstrat_creatinine, hr_creatinine, CAN_MOST_RECENT_CREAT),
    hr_albumin = coalesce(hr_riskstrat_albumin, hr_albumin, CAN_TOT_ALBUMIN),
    hr_bilirubin = coalesce(hr_riskstrat_bilirubin, hr_bilirubin, CAN_TOT_BILI),
    hr_sodium = coalesce(hr_riskstrat_sodium, hr_sodium, CAN_LAST_SERUM_SODIUM),
    hr_bnp_source = case_when(
      !is.na(hr_riskstrat_bnp) ~ "RiskStratDataHR",
      !is.na(hr_bnp) ~ "statjust_hr1a",
      TRUE ~ NA_character_
    ),
    hr_bnp_type_raw = if_else(!is.na(hr_riskstrat_bnp), hr_riskstrat_bnp_type, NA_character_),
    hr_bnp = coalesce(hr_riskstrat_bnp, hr_bnp),
    hr_bnp_type = case_when(
      str_detect(str_to_upper(str_squish(coalesce(hr_bnp_type_raw, ""))), "NT") ~ "NT-proBNP",
      str_detect(str_to_upper(str_squish(coalesce(hr_bnp_type_raw, ""))), "BNP") ~ "BNP",
      !is.na(hr_bnp) ~ "Unknown",
      TRUE ~ "Missing"
    ),
    hr_bnp_regular_indicator = as.integer(hr_bnp_type %in% c("BNP", "Unknown")),
    hr_bnp_ntpro_indicator = as.integer(hr_bnp_type == "NT-proBNP"),
    hr_albumin_observed = as.numeric(hr_albumin),
    hr_bilirubin_observed = as.numeric(hr_bilirubin),
    hr_sodium_observed = as.numeric(hr_sodium),
    hr_bnp_observed = as.numeric(hr_bnp),
    hr_short_mcs = coalesce(hr_short_mcs, flag_yes(CAN_ECMO), 0L),
    hr_durable_lvad = if_else(WL_ORG == "HR", flag_durable_vad_brand(CAN_VAD1, CAN_VAD2), 0L),
    lung_fev1_observed = as.numeric(CAN_FEV1),
    lung_fvc_observed = as.numeric(CAN_FVC),
    lung_pco2_observed = as.numeric(CAN_PCO2),
    lung_resting_o2_observed = as.numeric(CAN_AT_REST_O2),
    lung_six_min_walk_observed = as.numeric(CAN_SIX_MIN_WALK),
    lung_pulm_art_mean_observed = as.numeric(CAN_PULM_ART_MEAN),
    lung_cardiac_output_observed = as.numeric(CAN_CARDIAC_OUTPUT),
    lung_fev1_score_input = coalesce(lung_fev1_observed, 0),
    lung_fvc_score_input = coalesce(lung_fvc_observed, 0),
    lung_pco2_score_input = coalesce(lung_pco2_observed, 0),
    lung_resting_o2_score_input = coalesce(lung_resting_o2_observed, 0),
    lung_six_min_walk_score_input = coalesce(lung_six_min_walk_observed, 0),
    lung_pulm_art_mean_score_input = coalesce(lung_pulm_art_mean_observed, 0),
    lung_cardiac_output_score_input = coalesce(lung_cardiac_output_observed, 0),
    lung_ventilator_score_input = flag_yes(coalesce(as.character(CAN_VENTILATOR), as.character(CAN_ON_VENTILATOR))),
    lung_ecmo_score_input = flag_yes(CAN_ECMO),
    lung_corticosteroid_score_input = flag_yes(CAN_CORTICOST_DEPND)
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
    us_crs_proxy = 1.02 * hr_short_mcs +
      0.55 * safe_log1p(hr_bilirubin_imputed) -
      0.01 * hr_egfr +
      0.40 * safe_log1p(hr_bnp_imputed) * hr_bnp_regular_indicator +
      0.20 * safe_log1p(hr_bnp_imputed) * hr_bnp_ntpro_indicator -
      0.63 * hr_albumin_imputed -
      0.07 * hr_sodium_imputed -
      1.12 * hr_durable_lvad,
    lung_las_cas_component_score = rowSums(cbind(
      lung_fev1_score_input * -0.01,
      lung_fvc_score_input * -0.005,
      lung_pco2_score_input * 0.02,
      lung_resting_o2_score_input * 0.01,
      lung_six_min_walk_score_input * -0.002,
      lung_pulm_art_mean_score_input * 0.01,
      lung_cardiac_output_score_input * -0.05,
      lung_ventilator_score_input * 1.0,
      lung_ecmo_score_input * 1.5,
      lung_corticosteroid_score_input * 0.3
    ), na.rm = TRUE),
    organ_score = case_when(
      WL_ORG == "HR" ~ us_crs_proxy,
      WL_ORG == "KI" ~ NA_real_,
      WL_ORG == "LI" ~ liver_meld,
      WL_ORG == "LU" ~ lung_las_cas_component_score
    ),
    organ_score_name = case_when(
      WL_ORG == "HR" ~ "US-CRS proxy baseline",
      WL_ORG == "KI" ~ "Age + dialysis + diabetes baseline covariates",
      WL_ORG == "LI" ~ "Initial MELD, last if initial missing",
      WL_ORG == "LU" ~ "LAS/CAS component proxy baseline"
    )
  ) %>%
  filter(
    index_year >= analysis_start_year,
    index_year <= max(pollution_end_year),
    CAN_SOURCE %in% c("A", "R"),
    !is.na(index_date),
    !is.na(observed_end_date),
    observed_end_date >= index_date,
    !is.na(candidate_zip)
  )

registration_cohort <- cohort %>%
  arrange(PERS_ID, WL_ORG, index_date, observed_end_date, waitlist_row_id) %>%
  group_by(PERS_ID, WL_ORG) %>%
  mutate(person_organ_episode = 1L) %>%
  ungroup()

collapse_episode_audit <- registration_cohort %>%
  group_by(WL_ORG, PERS_ID, person_organ_episode) %>%
  summarise(
    episode_start = min(index_date, na.rm = TRUE),
    episode_end = max(observed_end_date, na.rm = TRUE),
    registrations = n(),
    centers = n_distinct(listing_center),
    multicenter = registrations > 1L & centers > 1L,
    .groups = "drop"
  )

collapse_summary <- collapse_episode_audit %>%
  group_by(WL_ORG) %>%
  summarise(
    registration_rows_before_collapse = sum(registrations),
    person_waitlist_episodes_after_collapse = n(),
    collapsed_episode_count = sum(registrations > 1L),
    collapsed_registration_rows = sum(registrations[registrations > 1L]),
    removed_duplicate_registration_rows = registration_rows_before_collapse - person_waitlist_episodes_after_collapse,
    multicenter_collapsed_episode_count = sum(multicenter),
    multicenter_collapsed_registration_rows = sum(registrations[multicenter]),
    .groups = "drop"
  ) %>%
  mutate(organ_label = recode(WL_ORG, !!!organ_labels), .before = 1)

write_csv(
  collapse_summary,
  file.path(out_dir, "overlapping_multilist_collapse_summary.csv")
)

log_msg(
  "Collapsed same-person/same-organ listings: ",
  sum(collapse_summary$registration_rows_before_collapse),
  " registrations -> ",
  sum(collapse_summary$person_waitlist_episodes_after_collapse),
  " person waitlist episodes; removed ",
  sum(collapse_summary$removed_duplicate_registration_rows),
  " duplicate registration rows"
)

registration_dt <- as.data.table(registration_cohort)
setorder(registration_dt, PERS_ID, WL_ORG, person_organ_episode, index_date, observed_end_date, waitlist_row_id)
episode_keys <- c("PERS_ID", "WL_ORG", "person_organ_episode")

baseline_dt <- registration_dt[, .SD[1], by = episode_keys]
episode_event_dt <- registration_dt[, {
  outcome_priority <- event_priority(adverse_event, transplant_or_improvement, other_exit)
  outcome_rows <- which(outcome_priority > 0L)
  if (length(outcome_rows) > 0L) {
    event_dat <- .SD[outcome_rows]
    event_dat[, outcome_priority := outcome_priority[outcome_rows]]
    setorder(event_dat, observed_end_date, outcome_priority, waitlist_row_id)
    final_event <- event_dat[.N]
    episode_end <- final_event$observed_end_date
    episode_adverse <- as.integer(final_event$adverse_event == 1L)
    episode_good <- as.integer(final_event$transplant_or_improvement == 1L)
    episode_other <- as.integer(final_event$other_exit == 1L)
  } else {
    episode_end <- max(observed_end_date, na.rm = TRUE)
    episode_adverse <- 0L
    episode_good <- 0L
    episode_other <- as.integer(any(other_exit == 1L, na.rm = TRUE))
  }
  episode_observed <- episode_adverse == 1L || episode_good == 1L || episode_other == 1L
  .(
    episode_observed_end_date = episode_end,
    episode_raw_event_date = if (episode_observed) episode_end else as.Date(NA),
    episode_event_observed = episode_observed,
    episode_adverse_event = episode_adverse,
    episode_transplant_or_improvement = episode_good,
    episode_other_exit = episode_other,
    episode_listing_count = .N,
    episode_center_count = uniqueN(listing_center),
    episode_px_ids = paste(unique(PX_ID), collapse = ";")
  )
}, by = episode_keys]

baseline_dt[
  episode_event_dt,
  on = episode_keys,
  `:=`(
    observed_end_date = i.episode_observed_end_date,
    raw_event_date = i.episode_raw_event_date,
    event_observed = i.episode_event_observed,
    adverse_event = i.episode_adverse_event,
    transplant_or_improvement = i.episode_transplant_or_improvement,
    other_exit = i.episode_other_exit,
    episode_listing_count = i.episode_listing_count,
    episode_center_count = i.episode_center_count,
    episode_px_ids = i.episode_px_ids
  )
]

cohort <- as_tibble(baseline_dt) %>%
  mutate(
    event_type = factor(
      case_when(
        adverse_event == 1L ~ "death_or_deterioration_delist",
        transplant_or_improvement == 1L ~ "transplant_or_improvement",
        other_exit == 1L ~ "other_exit",
        TRUE ~ "administrative_censor"
      ),
      levels = c("administrative_censor", "death_or_deterioration_delist", "transplant_or_improvement", "other_exit")
    )
  ) %>%
  arrange(WL_ORG, PERS_ID, index_date, observed_end_date, PX_ID) %>%
  mutate(
    source_waitlist_row_id = waitlist_row_id,
    waitlist_row_id = row_number(),
    index_year = as.integer(format(index_date, "%Y")),
    listing_year = factor(index_year)
  )

write_csv(
  cohort %>%
    select(
      waitlist_row_id, source_waitlist_row_id, PERS_ID, PX_ID, episode_px_ids,
      WL_ORG, person_organ_episode, episode_listing_count, episode_center_count,
      candidate_zip, index_date, observed_end_date,
      event_type, adverse_event, transplant_or_improvement, other_exit,
      age, sex, race, listing_year, organ_score, organ_score_name, listing_center,
      hr_riskstrat_creatinine, hr_riskstrat_albumin, hr_riskstrat_bilirubin,
      hr_riskstrat_sodium, hr_riskstrat_bnp, hr_bnp_source, hr_bnp_type_raw,
      hr_bnp_type, hr_bnp_regular_indicator, hr_bnp_ntpro_indicator,
      hr_albumin_observed, hr_bilirubin_observed, hr_sodium_observed, hr_bnp_observed,
      hr_albumin_imputed, hr_bilirubin_imputed, hr_egfr, hr_sodium_imputed,
      hr_durable_lvad, hr_short_mcs, hr_bnp_imputed,
      kidney_on_dialysis, kidney_no_dialysis_time, kidney_dialysis_years, kidney_diabetes, dialysis_years, liver_meld,
      lung_fev1_observed, lung_fvc_observed, lung_pco2_observed,
      lung_resting_o2_observed, lung_six_min_walk_observed,
      lung_pulm_art_mean_observed, lung_cardiac_output_observed,
      lung_fev1_score_input, lung_fvc_score_input, lung_pco2_score_input,
      lung_resting_o2_score_input, lung_six_min_walk_score_input,
      lung_pulm_art_mean_score_input, lung_cardiac_output_score_input,
      lung_ventilator_score_input, lung_ecmo_score_input, lung_corticosteroid_score_input
    ),
  file.path(out_dir, "primary_waitlist_period_pollution_analysis_dataset.csv.gz")
)

log_msg("Wrote the deduplicated cohort to ", out_dir)
