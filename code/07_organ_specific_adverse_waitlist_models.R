#!/usr/bin/env Rscript

runtime_source <- file.path("code", "r_runtime.R")
if (!file.exists(runtime_source)) runtime_source <- "r_runtime.R"
if (file.exists(runtime_source)) {
  source(runtime_source)
  ensure_user_library()
}
if (.Platform$OS.type == "windows") {
  r_minor <- strsplit(R.version$minor, "[.]")[[1]][[1]]
  user_lib <- Sys.getenv("R_LIBS_USER", unset = NA_character_)
  if (!is.na(user_lib) && nzchar(user_lib)) user_lib <- gsub("\\\\", "/", user_lib)
  if (is.na(user_lib) || !nzchar(user_lib)) {
    userprofile <- Sys.getenv("USERPROFILE", unset = path.expand("~"))
    user_lib <- paste0(gsub("\\\\", "/", userprofile), "/AppData/Local/R/win-library/", R.version$major, ".", r_minor)
  }
  .libPaths(c(user_lib, .libPaths()))
}
if (Sys.getenv("DEBUG_R_LIBPATHS", "FALSE") == "TRUE") {
  message("R_LIBS_USER=", Sys.getenv("R_LIBS_USER", unset = ""))
  message(".libPaths=", paste(.libPaths(), collapse = "; "))
}

suppressPackageStartupMessages({
  library(broom)
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(haven)
  library(readr)
  library(scales)
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
out_dir <- file.path("output", "organ_specific_adverse_waitlist")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
output_suffix <- Sys.getenv("OUTPUT_SUFFIX", unset = "")
if (nzchar(output_suffix) && !startsWith(output_suffix, "_")) output_suffix <- paste0("_", output_suffix)
step_dir <- file.path(out_dir, paste0("timevarying_steps", output_suffix))
dir.create(step_dir, recursive = TRUE, showWarnings = FALSE)
progress_path <- file.path(step_dir, "progress.csv")

analysis_start_year <- 2006L
analysis_end_year <- 2023L
target_organs <- c("HR", "KI", "LI", "LU")
target_organs_env <- Sys.getenv("TARGET_ORGANS", unset = "")
if (nzchar(target_organs_env)) {
  target_organs <- intersect(str_split(target_organs_env, ",", simplify = TRUE)[1, ], target_organs)
}

log_msg <- function(...) {
  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
  flush.console()
}

make_progress <- function(total, label) {
  if (total <= 0L) return(NULL)
  pb <- txtProgressBar(min = 0, max = total, style = 3)
  attr(pb, "label") <- label
  attr(pb, "total") <- total
  pb
}

tick_progress <- function(pb, value, detail = "") {
  if (is.null(pb)) return(invisible(NULL))
  setTxtProgressBar(pb, value)
  label <- attr(pb, "label")
  total <- attr(pb, "total")
  if (nzchar(detail)) message("\n", label, " ", value, "/", total, " | ", detail)
  invisible(NULL)
}

close_progress <- function(pb) {
  if (!is.null(pb)) close(pb)
}

safe_slug <- function(x) {
  str_to_lower(str_replace_all(x, "[^A-Za-z0-9]+", "_")) %>%
    str_replace_all("^_|_$", "")
}

append_progress <- function(row) {
  row <- as_tibble(row)
  progress_cols <- c(
    "organ", "model", "exposure_window", "exposure_form", "status", "error",
    "started_at", "finished_at", "elapsed_seconds", "n_intervals",
    "n_registrations", "adverse_events", "concordance"
  )
  for (col in setdiff(progress_cols, names(row))) row[[col]] <- NA
  row <- row[, progress_cols]
  write_csv(row, progress_path, append = file.exists(progress_path))
}

clean_zip <- function(x) {
  z <- str_extract(as.character(x), "[0-9]+")
  z <- if_else(nchar(z) > 5L, substr(z, 1, 5), z)
  if_else(is.na(z), NA_character_, str_pad(z, 5, side = "left", pad = "0"))
}

mode_value <- function(x) {
  ux <- na.omit(x)
  if (length(ux) == 0L) return(NA_character_)
  names(which.max(table(ux)))
}

first_existing <- function(data, vars, default = NA) {
  for (var in vars) {
    if (var %in% names(data)) return(data[[var]])
  }
  rep(default, nrow(data))
}

as_factor_missing <- function(x) {
  factor(if_else(is.na(x) | as.character(x) == "", "Missing", as.character(x)))
}

scale_01_to_pct <- function(x) {
  as.numeric(x) * 100
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
  "CAN_PRIMARY_PAY", "CAN_EDUCATION", "CAN_FUNCTN_STAT", "CAN_MED_COND",
  "CAN_LISTING_CTR_CD", "CAN_DGN", "CAN_DGN2",
  "CAN_INIT_STAT", "CAN_LAST_STAT", "CAN_INIT_ACT_STAT_CD",
  "CAN_LIFE_SUPPORT", "CAN_ECMO", "CAN_VENTILATOR", "CAN_IABP", "CAN_IV_INOTROP",
  "CAN_VAD_TAH", "CAN_VAD1", "CAN_VAD2", "CAN_TAH", "CAN_VASC_ASSIST", "CAN_ICU",
  "CAN_INOTROP", "CAN_ON_VENTILATOR",
  "CAN_MOST_RECENT_CREAT", "CAN_TOT_ALBUMIN", "CAN_DIAB", "CAN_DIAB_TY",
  "CAN_ON_DIAL", "CAN_DIAL", "CAN_DIAL_DT", "CAN_GFR", "CAN_INIT_ALLOC_PRA",
  "CAN_INIT_CUR_PRA", "CAN_INIT_SRTR_PEAK_PRA", "CAN_LATEST_SRTR_PEAK_PRA",
  "CAN_INIT_SRTR_LAB_MELD", "CAN_INIT_SRTR_LAB_MELD_TY", "CAN_LAST_SRTR_LAB_MELD",
  "CAN_TOT_BILI", "CAN_CTP_SCORE", "CAN_LAST_INR", "CAN_LAST_SERUM_SODIUM",
  "CAN_LAST_DIAL_PRIOR_WEEK", "CAN_ASCITES", "CAN_ENCEPH",
  "CAN_FVC", "CAN_FEV1", "CAN_PCO2", "CAN_AT_REST_O2", "CAN_SIX_MIN_WALK",
  "CAN_PULM_ART_MEAN", "CAN_PCW_MEAN", "CAN_CARDIAC_OUTPUT", "CAN_CORTICOST_DEPND"
)

log_msg("Reading candidate ZIP supplement")
candidate_zip <- read_sas(file.path(supp_dir, "canzip2506.sas7bdat")) %>%
  transmute(PERS_ID, PX_ID, candidate_zip = clean_zip(CAN_PERM_ZIP)) %>%
  distinct(PERS_ID, PX_ID, .keep_all = TRUE)

log_msg("Reading pollution exposure tables")
pollution_annual <- read_pollution_annual() %>%
  filter(year >= analysis_start_year - 1L, year <= analysis_end_year)

community_path <- file.path(community_dir, "zcta_acs_community_covariates_2005_2023.csv.gz")
if (file.exists(community_path)) {
  log_msg("Reading community covariates")
  community <- read_csv(community_path, col_types = cols(zip = col_character(), .default = col_guess()), show_col_types = FALSE) %>%
    transmute(
      zip = clean_zip(zip),
      analysis_year,
      median_household_income,
      pct_poverty = scale_01_to_pct(pct_poverty),
      pct_bachelor_plus = scale_01_to_pct(pct_bachelor_plus),
      pct_unemployed = scale_01_to_pct(pct_unemployed),
      pct_no_vehicle = scale_01_to_pct(pct_no_vehicle),
      pct_nonwhite = scale_01_to_pct(pct_nonwhite),
      pct_black = scale_01_to_pct(pct_black),
      pct_hispanic = scale_01_to_pct(pct_hispanic)
    )
} else {
  log_msg("Community covariates not found; run code/06_build_community_covariates.R. Models will omit community adjustment.")
  community <- NULL
}

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

log_msg("Constructing heart/kidney/liver/lung cohort")
cohort <- candidate_data %>%
  filter(WL_ORG %in% target_organs) %>%
  left_join(candidate_zip, by = c("PERS_ID", "PX_ID")) %>%
  mutate(
    waitlist_row_id = row_number(),
    index_date = coalesce(CAN_ACTIVATE_DT, CAN_LISTING_DT),
    index_year = as.integer(format(index_date, "%Y")),
    adverse_waitlist_event = CAN_REM_CD %in% c(8, 13),
    death_event = coalesce(CAN_REM_CD == 8, FALSE),
    deterioration_event = coalesce(CAN_REM_CD == 13, FALSE),
    improved_delist = coalesce(CAN_REM_CD == 12, FALSE),
    medically_unsuitable = coalesce(CAN_REM_CD == 5, FALSE),
    transplant_event = CAN_REM_CD %in% c(4, 14, 15, 18, 19, 21, 22, 23),
    other_delist_event = !adverse_waitlist_event & !transplant_event & !improved_delist & !is.na(CAN_REM_CD),
    event_date = coalesce(CAN_REM_DT, CAN_DEATH_DT, CAN_ENDWLFU),
    end_date = event_date,
    followup_days = as.numeric(end_date - index_date),
    age = CAN_AGE_AT_LISTING,
    sex = as_factor_missing(if_else(CAN_GENDER %in% c("M", "F"), CAN_GENDER, NA_character_)),
    race_srtr = as_factor_missing(CAN_RACE_SRTR),
    ethnicity_srtr = as_factor_missing(CAN_ETHNICITY_SRTR),
    abo = as_factor_missing(CAN_ABO),
    primary_pay = as_factor_missing(CAN_PRIMARY_PAY),
    education = as_factor_missing(CAN_EDUCATION),
    functional_status = as_factor_missing(CAN_FUNCTN_STAT),
    medical_condition = as_factor_missing(CAN_MED_COND),
    diagnosis = as_factor_missing(CAN_DGN),
    init_status = as_factor_missing(CAN_INIT_STAT),
    listing_center = as_factor_missing(CAN_LISTING_CTR_CD),
    index_year_factor = factor(index_year),
    bmi = CAN_BMI,
    creatinine = CAN_MOST_RECENT_CREAT,
    albumin = CAN_TOT_ALBUMIN,
    diabetes = as_factor_missing(CAN_DIAB),
    life_support = as_factor_missing(CAN_LIFE_SUPPORT),
    ecmo = as_factor_missing(CAN_ECMO),
    ventilator = as_factor_missing(coalesce(as.character(CAN_VENTILATOR), as.character(CAN_ON_VENTILATOR))),
    iabp = as_factor_missing(CAN_IABP),
    inotropes = as_factor_missing(coalesce(as.character(CAN_IV_INOTROP), as.character(CAN_INOTROP))),
    vad_tah = as_factor_missing(CAN_VAD_TAH),
    vad1 = as_factor_missing(CAN_VAD1),
    vad2 = as_factor_missing(CAN_VAD2),
    icu = as_factor_missing(CAN_ICU),
    on_dialysis = as_factor_missing(coalesce(as.character(CAN_ON_DIAL), as.character(CAN_DIAL))),
    dialysis_years = as.numeric(index_date - CAN_DIAL_DT) / 365.25,
    dialysis_years = if_else(is.na(dialysis_years) | dialysis_years < 0, NA_real_, dialysis_years),
    dialysis_years_zero = coalesce(dialysis_years, 0),
    log1p_dialysis_years = log1p(dialysis_years_zero),
    gfr = CAN_GFR,
    alloc_pra = CAN_INIT_ALLOC_PRA,
    current_pra = CAN_INIT_CUR_PRA,
    peak_pra = coalesce(CAN_INIT_SRTR_PEAK_PRA, CAN_LATEST_SRTR_PEAK_PRA),
    init_lab_meld = CAN_INIT_SRTR_LAB_MELD,
    init_meld_type = as_factor_missing(CAN_INIT_SRTR_LAB_MELD_TY),
    bilirubin = CAN_TOT_BILI,
    ctp_score = CAN_CTP_SCORE,
    inr = CAN_LAST_INR,
    sodium = CAN_LAST_SERUM_SODIUM,
    ascites = as_factor_missing(CAN_ASCITES),
    encephalopathy = as_factor_missing(CAN_ENCEPH),
    fvc = CAN_FVC,
    fev1 = CAN_FEV1,
    pco2 = CAN_PCO2,
    rest_o2 = CAN_AT_REST_O2,
    six_min_walk = CAN_SIX_MIN_WALK,
    pulm_art_mean = CAN_PULM_ART_MEAN,
    wedge_mean = CAN_PCW_MEAN,
    cardiac_output = CAN_CARDIAC_OUTPUT,
    corticosteroid_dependent = as_factor_missing(CAN_CORTICOST_DEPND)
  ) %>%
  filter(
    index_year >= analysis_start_year,
    index_year <= analysis_end_year,
    CAN_SOURCE %in% c("A", "R"),
    !is.na(index_date),
    !is.na(end_date),
    !is.na(followup_days),
    followup_days >= 0,
    !is.na(candidate_zip)
  ) %>%
  mutate(
    followup_days = pmax(followup_days, 0.5),
    adverse_waitlist_event = as.integer(adverse_waitlist_event),
    death_event = as.integer(death_event),
    deterioration_event = as.integer(deterioration_event),
    transplant_event = as.integer(transplant_event),
    improved_delist = as.integer(improved_delist),
    other_delist_event = as.integer(other_delist_event)
  )

baseline_exposure <- pollution_annual %>%
  transmute(
    candidate_zip = zip,
    prior_exposure_year = year,
    pm25_prior_1y_5ug = pm25_5ug,
    o3_prior_1y_10ppb = o3_10ppb,
    no2_prior_1y_10unit = no2_10unit
  )

cohort <- cohort %>%
  mutate(prior_exposure_year = index_year - 1L) %>%
  left_join(baseline_exposure, by = c("candidate_zip", "prior_exposure_year"))

if (!is.null(community)) {
  cohort <- cohort %>%
    left_join(community, by = c("candidate_zip" = "zip", "index_year" = "analysis_year"))
}

base_required <- c(
  "pm25_prior_1y_5ug", "o3_prior_1y_10ppb", "no2_prior_1y_10unit",
  "age", "sex", "race_srtr", "ethnicity_srtr", "bmi", "abo", "listing_center",
  "index_year_factor"
)

write_csv(
  cohort %>%
    group_by(WL_ORG) %>%
    summarise(
      n = n(),
      people = n_distinct(PERS_ID),
      adverse_events = sum(adverse_waitlist_event),
      deaths = sum(death_event),
      deteriorations = sum(deterioration_event),
      transplants = sum(transplant_event),
      improved_delist = sum(improved_delist),
      medically_unsuitable = sum(medically_unsuitable, na.rm = TRUE),
      complete_baseline_exposure = sum(complete.cases(across(all_of(base_required)))),
      .groups = "drop"
    ),
  file.path(out_dir, "revised_cohort_summary_by_organ.csv")
)

community_terms <- if (!is.null(community)) {
  c(
    "ns(median_household_income, df = 3)",
    "ns(pct_poverty, df = 3)",
    "ns(pct_bachelor_plus, df = 3)",
    "ns(pct_unemployed, df = 3)",
    "ns(pct_no_vehicle, df = 3)",
    "ns(pct_black, df = 3)",
    "ns(pct_hispanic, df = 3)"
  )
} else {
  character()
}

organ_terms <- list(
  HR = c(
    "init_status", "medical_condition", "diagnosis", "functional_status",
    "life_support", "ecmo", "ventilator", "iabp", "inotropes", "vad_tah",
    "vad1", "vad2", "ns(creatinine, df = 3)", "ns(albumin, df = 3)", "diabetes"
  ),
  KI = c(
    "medical_condition", "diagnosis", "functional_status", "on_dialysis",
    "log1p_dialysis_years", "ns(gfr, df = 3)", "ns(creatinine, df = 3)",
    "alloc_pra", "current_pra", "peak_pra", "diabetes"
  ),
  LI = c(
    "init_status", "medical_condition", "diagnosis", "functional_status",
    "ns(init_lab_meld, df = 4)", "init_meld_type", "ns(bilirubin, df = 3)",
    "ns(creatinine, df = 3)", "ns(albumin, df = 3)", "ns(inr, df = 3)",
    "ns(sodium, df = 3)", "ascites", "encephalopathy", "diabetes"
  ),
  LU = c(
    "init_status", "medical_condition", "diagnosis", "functional_status",
    "life_support", "ecmo", "ventilator", "icu", "ns(fvc, df = 3)",
    "ns(fev1, df = 3)", "ns(pco2, df = 3)", "ns(rest_o2, df = 3)",
    "ns(six_min_walk, df = 3)", "ns(pulm_art_mean, df = 3)",
    "ns(wedge_mean, df = 3)", "ns(cardiac_output, df = 3)",
    "corticosteroid_dependent", "diabetes"
  )
)

common_terms <- c(
  "pm25_prior_1y_5ug", "o3_prior_1y_10ppb", "no2_prior_1y_10unit",
  "ns(age, df = 4)", "sex", "race_srtr", "ethnicity_srtr", "ns(bmi, df = 3)",
  "abo", "primary_pay", "education",
  community_terms,
  "strata(listing_center)", "strata(index_year_factor)", "cluster(PERS_ID)"
)

extract_formula_vars <- function(terms) {
  vars <- unique(unlist(str_extract_all(terms, "[A-Za-z][A-Za-z0-9_]*")))
  setdiff(vars, c("ns", "df", "strata", "cluster"))
}

term_is_usable <- function(term, dat, min_observed = 0.60) {
  vars <- intersect(extract_formula_vars(term), names(dat))
  if (length(vars) == 0L) return(TRUE)

  for (var in vars) {
    x <- dat[[var]]
    observed <- mean(!is.na(x))
    if (observed < min_observed && !is.factor(x)) return(FALSE)
    ux <- unique(na.omit(x))
    if (length(ux) < 2L && !grepl("^cluster\\(", term)) return(FALSE)
    if (grepl("^ns\\(", term) && length(ux) < 5L) return(FALSE)
  }

  TRUE
}

filter_usable_terms <- function(terms, dat, protected_terms = character()) {
  keep <- vapply(terms, function(term) term %in% protected_terms || term_is_usable(term, dat), logical(1))
  dropped <- terms[!keep]
  if (length(dropped) > 0L) {
    log_msg("Dropping low-information terms: ", paste(dropped, collapse = ", "))
  }
  terms[keep]
}

fit_organ_model <- function(org, dat, extra_terms, exposure_terms, model_name) {
  rhs <- c(exposure_terms, setdiff(c(common_terms, extra_terms), c("pm25_prior_1y_5ug", "o3_prior_1y_10ppb", "no2_prior_1y_10unit")))
  rhs <- filter_usable_terms(rhs, dat, protected_terms = exposure_terms)
  vars_needed <- intersect(extract_formula_vars(rhs), names(dat))
  model_dat <- dat %>%
    filter(complete.cases(across(all_of(vars_needed)))) %>%
    droplevels()

  if (nrow(model_dat) == 0L || sum(model_dat$adverse_waitlist_event) < 100L) return(NULL)

  form <- as.formula(paste("Surv(followup_days, adverse_waitlist_event) ~", paste(rhs, collapse = " + ")))
  log_msg("Fitting ", model_name, " for ", org, " n=", nrow(model_dat), " adverse=", sum(model_dat$adverse_waitlist_event))
  fit <- coxph(form, data = model_dat, ties = "efron", robust = TRUE, x = FALSE, y = FALSE)

  tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
    mutate(
      organ = org,
      model = model_name,
      n = nrow(model_dat),
      adverse_events = sum(model_dat$adverse_waitlist_event),
      concordance = unname(summary(fit)$concordance[1]),
      .before = 1
    )
}

baseline_results <- list()
pollution_results <- list()

if (Sys.getenv("RUN_BASELINE_MODELS", "TRUE") != "FALSE") {
  for (org in target_organs) {
    dat <- cohort %>% filter(WL_ORG == org)
    extra <- organ_terms[[org]]
    baseline_results[[org]] <- fit_organ_model(org, dat, extra, character(), "clinical_center_community")
    pollution_results[[org]] <- fit_organ_model(
      org,
      dat,
      extra,
      c("pm25_prior_1y_5ug", "o3_prior_1y_10ppb", "no2_prior_1y_10unit"),
      "clinical_center_community_plus_prior_pollution"
    )
  }

  bind_rows(c(baseline_results, pollution_results)) %>%
    write_csv(file.path(out_dir, "organ_specific_baseline_prior_pollution_cox_results.csv"))
} else {
  log_msg("Skipping baseline organ models because RUN_BASELINE_MODELS=FALSE")
}

log_msg("Constructing lagged annual waitlist exposure intervals")
intervals <- as.data.table(cohort %>%
  filter(!is.na(pm25_prior_1y_5ug), !is.na(o3_prior_1y_10ppb), !is.na(no2_prior_1y_10unit)) %>%
  select(
    waitlist_row_id, PERS_ID, PX_ID, WL_ORG, candidate_zip, index_date, end_date, followup_days,
    adverse_waitlist_event, death_event, deterioration_event, transplant_event, improved_delist,
    age, sex, race_srtr, ethnicity_srtr, bmi, abo, primary_pay, education, listing_center,
    index_year, index_year_factor, all_of(names(cohort)[names(cohort) %in% extract_formula_vars(unlist(organ_terms))]),
    pm25_prior_1y_5ug, o3_prior_1y_10ppb, no2_prior_1y_10unit,
    any_of(c(
      "median_household_income", "pct_poverty", "pct_bachelor_plus", "pct_unemployed",
      "pct_no_vehicle", "pct_black", "pct_hispanic"
    ))
  ))

intervals[, interval_years := lapply(seq_len(.N), function(i) {
  start_year <- as.integer(format(index_date[i], "%Y"))
  stop_year <- min(as.integer(format(end_date[i], "%Y")), analysis_end_year)
  if (is.na(start_year) || is.na(stop_year) || stop_year < start_year) return(integer())
  seq.int(start_year, stop_year)
})]
intervals <- intervals[lengths(interval_years) > 0L]
interval_lengths <- lengths(intervals$interval_years)
interval_row_index <- rep(seq_len(nrow(intervals)), interval_lengths)
interval_year_values <- unlist(intervals$interval_years, use.names = FALSE)
intervals <- intervals[interval_row_index]
intervals[, interval_year := interval_year_values]
intervals[, interval_years := NULL]
intervals[, interval_start_date := pmax(index_date, as.Date(paste0(interval_year, "-01-01")))]
intervals[, interval_end_date := pmin(end_date, as.Date(paste0(interval_year, "-12-31")))]
intervals <- intervals[interval_end_date >= interval_start_date]
intervals[, tstart := as.numeric(interval_start_date - index_date)]
intervals[, tstop := pmax(as.numeric(interval_end_date - index_date), tstart + 0.5)]
intervals[, adverse_event_interval := as.integer(adverse_waitlist_event == 1L & interval_end_date >= end_date)]
intervals[, exposure_year := interval_year - 1L]

tv_exp <- as.data.table(pollution_annual %>%
  transmute(
    candidate_zip = zip,
    exposure_year = year,
    pm25_waitlist_5ug = pm25_5ug,
    o3_waitlist_10ppb = o3_10ppb,
    no2_waitlist_10unit = no2_10unit
  ))
setkey(intervals, candidate_zip, exposure_year)
setkey(tv_exp, candidate_zip, exposure_year)
intervals <- tv_exp[intervals]
intervals <- intervals[!is.na(pm25_waitlist_5ug) & !is.na(o3_waitlist_10ppb) & !is.na(no2_waitlist_10unit)]
intervals[, interval_days := pmax(tstop - tstart, 0.5)]
setorder(intervals, waitlist_row_id, tstart, tstop)
intervals[, pm25_waitlist_cum_5ug := shift(cumsum(pm25_waitlist_5ug * interval_days), fill = NA_real_) / shift(cumsum(interval_days), fill = NA_real_), by = waitlist_row_id]
intervals[, o3_waitlist_cum_10ppb := shift(cumsum(o3_waitlist_10ppb * interval_days), fill = NA_real_) / shift(cumsum(interval_days), fill = NA_real_), by = waitlist_row_id]
intervals[, no2_waitlist_cum_10unit := shift(cumsum(no2_waitlist_10unit * interval_days), fill = NA_real_) / shift(cumsum(interval_days), fill = NA_real_), by = waitlist_row_id]
intervals[is.na(pm25_waitlist_cum_5ug), pm25_waitlist_cum_5ug := pm25_prior_1y_5ug]
intervals[is.na(o3_waitlist_cum_10ppb), o3_waitlist_cum_10ppb := o3_prior_1y_10ppb]
intervals[is.na(no2_waitlist_cum_10unit), no2_waitlist_cum_10unit := no2_prior_1y_10unit]

write_csv(
  as_tibble(intervals) %>%
    select(
      waitlist_row_id, PERS_ID, PX_ID, WL_ORG, candidate_zip, tstart, tstop,
      interval_year, exposure_year, adverse_event_interval, pm25_waitlist_5ug,
      o3_waitlist_10ppb, no2_waitlist_10unit, pm25_waitlist_cum_5ug,
      o3_waitlist_cum_10ppb, no2_waitlist_cum_10unit
    ) %>%
    slice_head(n = 10000),
  file.path(out_dir, "lagged_waitlist_exposure_intervals_qc_sample.csv")
)

write_csv(
  as_tibble(intervals) %>%
    group_by(WL_ORG) %>%
    summarise(
      n_intervals = n(),
      n_registrations = n_distinct(waitlist_row_id),
      adverse_events = sum(adverse_event_interval),
      min_exposure_year = min(exposure_year),
      max_exposure_year = max(exposure_year),
      .groups = "drop"
    ),
  file.path(out_dir, "lagged_waitlist_exposure_interval_summary.csv")
)

pollutant_definitions <- tribble(
  ~pollutant, ~current_term, ~cumulative_term,
  "PM2.5", "pm25_waitlist_5ug", "pm25_waitlist_cum_5ug",
  "Ozone", "o3_waitlist_10ppb", "o3_waitlist_cum_10ppb",
  "NO2", "no2_waitlist_10unit", "no2_waitlist_cum_10unit"
)

wald_test_terms <- function(fit, pattern) {
  coef_names <- names(coef(fit))
  idx <- grep(pattern, coef_names, fixed = FALSE)
  if (length(idx) == 0L) return(tibble(wald_chisq = NA_real_, df = NA_integer_, p_value = NA_real_))

  beta <- coef(fit)[idx]
  v <- vcov(fit)[idx, idx, drop = FALSE]
  stat <- as.numeric(t(beta) %*% solve(v, beta))
  tibble(wald_chisq = stat, df = length(idx), p_value = pchisq(stat, df = length(idx), lower.tail = FALSE))
}

contrast_iqr <- function(fit, model_dat, exposure_term, pollutant, organ, model_name, exposure_window, exposure_form) {
  exposure_value <- model_dat[[exposure_term]]
  q <- quantile(exposure_value, probs = c(0.25, 0.75), na.rm = TRUE)
  if (!all(is.finite(q)) || q[[1]] == q[[2]]) return(tibble())

  base_row <- model_dat[which.min(abs(exposure_value - median(exposure_value, na.rm = TRUE))), , drop = FALSE]
  low_row <- base_row
  high_row <- base_row
  low_row[[exposure_term]] <- q[[1]]
  high_row[[exposure_term]] <- q[[2]]

  pred <- predict(fit, newdata = bind_rows(low_row, high_row), type = "lp", se.fit = TRUE)
  log_hr <- pred$fit[[2]] - pred$fit[[1]]
  se <- sqrt(pred$se.fit[[1]]^2 + pred$se.fit[[2]]^2)

  pattern <- paste0("(^|\\`|ns\\()", exposure_term)
  wald <- wald_test_terms(fit, pattern)

  tibble(
    organ = organ,
    pollutant = pollutant,
    model = model_name,
    exposure_window = exposure_window,
    exposure_form = exposure_form,
    contrast = "p75_vs_p25",
    p25 = as.numeric(q[[1]]),
    p75 = as.numeric(q[[2]]),
    hazard_ratio = exp(log_hr),
    conf_low = exp(log_hr - 1.96 * se),
    conf_high = exp(log_hr + 1.96 * se),
    n_intervals = nrow(model_dat),
    n_registrations = n_distinct(model_dat$waitlist_row_id),
    adverse_events = sum(model_dat$adverse_event_interval),
    concordance = unname(summary(fit)$concordance[1])
  ) %>%
    bind_cols(wald)
}

make_spline_curve <- function(fit, model_dat, exposure_term, pollutant, organ, model_name, exposure_window, exposure_form) {
  exposure_value <- model_dat[[exposure_term]]
  qs <- quantile(exposure_value, probs = c(0.01, 0.50, 0.99), na.rm = TRUE)
  if (!all(is.finite(qs)) || qs[[1]] == qs[[3]]) return(tibble())

  grid <- seq(qs[[1]], qs[[3]], length.out = 80)
  ref_row <- model_dat[which.min(abs(exposure_value - qs[[2]])), , drop = FALSE]
  newdat <- ref_row[rep(1, length(grid)), , drop = FALSE]
  newdat[[exposure_term]] <- grid

  ref_dat <- ref_row
  ref_dat[[exposure_term]] <- qs[[2]]
  pred <- predict(fit, newdata = bind_rows(ref_dat, newdat), type = "lp", se.fit = TRUE)
  ref_fit <- pred$fit[[1]]
  curve_fit <- pred$fit[-1]
  curve_se <- pred$se.fit[-1]

  tibble(
    organ = organ,
    pollutant = pollutant,
    model = model_name,
    exposure_window = exposure_window,
    exposure_form = exposure_form,
    exposure_value = grid,
    reference_value = as.numeric(qs[[2]]),
    hazard_ratio = exp(curve_fit - ref_fit),
    conf_low = exp(curve_fit - ref_fit - 1.96 * curve_se),
    conf_high = exp(curve_fit - ref_fit + 1.96 * curve_se),
    n_intervals = nrow(model_dat),
    n_registrations = n_distinct(model_dat$waitlist_row_id),
    adverse_events = sum(model_dat$adverse_event_interval),
    concordance = unname(summary(fit)$concordance[1])
  )
}

fit_timevarying_spec <- function(org, dat, extra, exposure_terms, protected_terms, model_name, exposure_window, exposure_form, pollutants) {
  step_started <- Sys.time()
  rhs <- c(exposure_terms, setdiff(c(common_terms, extra), c("pm25_prior_1y_5ug", "o3_prior_1y_10ppb", "no2_prior_1y_10unit")))
  rhs <- filter_usable_terms(rhs, dat, protected_terms = protected_terms)
  vars_needed <- intersect(extract_formula_vars(rhs), names(dat))
  model_dat <- dat %>%
    filter(complete.cases(across(all_of(vars_needed)))) %>%
    droplevels()
  if (nrow(model_dat) == 0L || sum(model_dat$adverse_event_interval) < 100L) {
    append_progress(tibble(
      organ = org,
      model = model_name,
      exposure_window = exposure_window,
      exposure_form = exposure_form,
      status = "skipped",
      started_at = as.character(step_started),
      finished_at = as.character(Sys.time()),
      elapsed_seconds = as.numeric(difftime(Sys.time(), step_started, units = "secs")),
      n_intervals = nrow(model_dat),
      adverse_events = sum(model_dat$adverse_event_interval)
    ))
    return(NULL)
  }

  form <- as.formula(paste("Surv(tstart, tstop, adverse_event_interval) ~", paste(rhs, collapse = " + ")))
  log_msg("START ", org, " | ", model_name, " | intervals=", nrow(model_dat), " | adverse=", sum(model_dat$adverse_event_interval))
  fit <- tryCatch(
    coxph(form, data = model_dat, ties = "efron", robust = TRUE, x = FALSE, y = FALSE),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    log_msg("ERROR ", org, " | ", model_name, " | ", conditionMessage(fit))
    append_progress(tibble(
      organ = org,
      model = model_name,
      exposure_window = exposure_window,
      exposure_form = exposure_form,
      status = "error",
      error = conditionMessage(fit),
      started_at = as.character(step_started),
      finished_at = as.character(Sys.time()),
      elapsed_seconds = as.numeric(difftime(Sys.time(), step_started, units = "secs")),
      n_intervals = nrow(model_dat),
      adverse_events = sum(model_dat$adverse_event_interval)
    ))
    return(NULL)
  }

  linear_results <- tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(term %in% protected_terms) %>%
    mutate(
      organ = org,
      model = model_name,
      exposure_window = exposure_window,
      exposure_form = exposure_form,
      pollutant = recode(term, !!!setNames(pollutants$pollutant, pollutants$term)),
      n_intervals = nrow(model_dat),
      n_registrations = n_distinct(model_dat$waitlist_row_id),
      adverse_events = sum(model_dat$adverse_event_interval),
      concordance = unname(summary(fit)$concordance[1]),
      .before = 1
    )

  contrast_results <- bind_rows(lapply(seq_len(nrow(pollutants)), function(i) {
    contrast_iqr(
      fit = fit,
      model_dat = model_dat,
      exposure_term = pollutants$term[[i]],
      pollutant = pollutants$pollutant[[i]],
      organ = org,
      model_name = model_name,
      exposure_window = exposure_window,
      exposure_form = exposure_form
    )
  }))

  curve_results <- if (exposure_form == "spline") {
    bind_rows(lapply(seq_len(nrow(pollutants)), function(i) {
      make_spline_curve(
        fit = fit,
        model_dat = model_dat,
        exposure_term = pollutants$term[[i]],
        pollutant = pollutants$pollutant[[i]],
        organ = org,
        model_name = model_name,
        exposure_window = exposure_window,
        exposure_form = exposure_form
      )
    }))
  } else {
    tibble()
  }

  step_slug <- paste(safe_slug(org), safe_slug(model_name), sep = "_")
  write_csv(linear_results, file.path(step_dir, paste0(step_slug, "_linear.csv")))
  write_csv(contrast_results, file.path(step_dir, paste0(step_slug, "_contrast.csv")))
  if (nrow(curve_results) > 0L) {
    write_csv(curve_results, file.path(step_dir, paste0(step_slug, "_curve.csv")))
  }
  elapsed <- as.numeric(difftime(Sys.time(), step_started, units = "secs"))
  log_msg("DONE  ", org, " | ", model_name, " | elapsed=", round(elapsed, 1), " sec")
  append_progress(tibble(
    organ = org,
    model = model_name,
    exposure_window = exposure_window,
    exposure_form = exposure_form,
    status = "done",
    started_at = as.character(step_started),
    finished_at = as.character(Sys.time()),
    elapsed_seconds = elapsed,
    n_intervals = nrow(model_dat),
    n_registrations = n_distinct(model_dat$waitlist_row_id),
    adverse_events = sum(model_dat$adverse_event_interval),
    concordance = unname(summary(fit)$concordance[1])
  ))

  list(linear = linear_results, contrast = contrast_results, curve = curve_results)
}

fit_organ_tv_model <- function(org) {
  dat <- as_tibble(intervals[WL_ORG == org])
  extra <- organ_terms[[org]]

  specs <- list()
  tv_windows <- str_split(Sys.getenv("TV_EXPOSURE_WINDOWS", "cumulative,lagged_annual"), ",", simplify = TRUE)[1, ]
  tv_windows <- intersect(tv_windows, c("cumulative", "lagged_annual"))
  tv_forms <- str_split(Sys.getenv("TV_EXPOSURE_FORMS", "linear,spline"), ",", simplify = TRUE)[1, ]
  tv_forms <- intersect(tv_forms, c("linear", "spline"))
  tv_model_types <- str_split(Sys.getenv("TV_MODEL_TYPES", "single,combined"), ",", simplify = TRUE)[1, ]
  tv_model_types <- intersect(tv_model_types, c("single", "combined"))

  for (window in tv_windows) {
    term_col <- if (window == "cumulative") "cumulative_term" else "current_term"
    window_terms <- pollutant_definitions[[term_col]]

    if ("single" %in% tv_model_types) {
      for (i in seq_len(nrow(pollutant_definitions))) {
        pollutant_name <- pollutant_definitions$pollutant[[i]]
        term <- window_terms[[i]]
        if ("linear" %in% tv_forms) {
          specs[[length(specs) + 1L]] <- list(
            model_name = paste0("single_pollutant_", str_to_lower(str_replace_all(pollutant_name, "[^A-Za-z0-9]+", "")), "_", window, "_linear"),
            exposure_window = window,
            exposure_form = "linear",
            exposure_terms = term,
            protected_terms = term,
            pollutants = tibble(pollutant = pollutant_name, term = term)
          )
        }
        if ("spline" %in% tv_forms) {
          specs[[length(specs) + 1L]] <- list(
            model_name = paste0("single_pollutant_", str_to_lower(str_replace_all(pollutant_name, "[^A-Za-z0-9]+", "")), "_", window, "_spline"),
            exposure_window = window,
            exposure_form = "spline",
            exposure_terms = paste0("ns(", term, ", df = 3)"),
            protected_terms = term,
            pollutants = tibble(pollutant = pollutant_name, term = term)
          )
        }
      }
    }

    if ("combined" %in% tv_model_types && "linear" %in% tv_forms) {
      specs[[length(specs) + 1L]] <- list(
        model_name = paste0("combined_pollutants_", window, "_linear"),
        exposure_window = window,
        exposure_form = "linear",
        exposure_terms = window_terms,
        protected_terms = window_terms,
        pollutants = tibble(pollutant = pollutant_definitions$pollutant, term = window_terms)
      )
    }
    if ("combined" %in% tv_model_types && "spline" %in% tv_forms) {
      specs[[length(specs) + 1L]] <- list(
        model_name = paste0("combined_pollutants_", window, "_spline"),
        exposure_window = window,
        exposure_form = "spline",
        exposure_terms = paste0("ns(", window_terms, ", df = 3)"),
        protected_terms = window_terms,
        pollutants = tibble(pollutant = pollutant_definitions$pollutant, term = window_terms)
      )
    }
  }

  pb <- make_progress(length(specs), paste0("Time-varying model grid for ", org))
  on.exit(close_progress(pb), add = TRUE)
  lapply(seq_along(specs), function(i) {
    spec <- specs[[i]]
    tick_progress(pb, i, spec$model_name)
    fit_timevarying_spec(
      org = org,
      dat = dat,
      extra = extra,
      exposure_terms = spec$exposure_terms,
      protected_terms = spec$protected_terms,
      model_name = spec$model_name,
      exposure_window = spec$exposure_window,
      exposure_form = spec$exposure_form,
      pollutants = spec$pollutants
    )
  })
}

if (Sys.getenv("RUN_TIMEVARYING_MODELS", "TRUE") != "FALSE") {
  organ_pb <- make_progress(length(target_organs), "Time-varying organ progress")
  tv_results <- unlist(lapply(seq_along(target_organs), function(i) {
    org <- target_organs[[i]]
    tick_progress(organ_pb, i, org)
    fit_organ_tv_model(org)
  }), recursive = FALSE)
  close_progress(organ_pb)
  tv_results <- Filter(Negate(is.null), tv_results)
  bind_rows(lapply(tv_results, function(x) x[["linear"]])) %>%
    write_csv(file.path(out_dir, paste0("organ_specific_timevarying_pollution_linear_cox_results", output_suffix, ".csv")))
  bind_rows(lapply(tv_results, function(x) x[["contrast"]])) %>%
    write_csv(file.path(out_dir, paste0("organ_specific_timevarying_pollution_spline_iqr_contrasts", output_suffix, ".csv")))
  bind_rows(lapply(tv_results, function(x) x[["curve"]])) %>%
    write_csv(file.path(out_dir, paste0("organ_specific_timevarying_pollution_spline_curves", output_suffix, ".csv")))
} else {
  log_msg("Skipping time-varying Cox models because RUN_TIMEVARYING_MODELS=FALSE")
}

log_msg("Wrote revised organ-specific adverse waitlist outputs to ", out_dir)
