#!/usr/bin/env Rscript

runtime_source <- file.path("code", "r_runtime.R")
if (!file.exists(runtime_source)) runtime_source <- "r_runtime.R"
if (file.exists(runtime_source)) {
  source(runtime_source)
  ensure_user_library()
}

suppressPackageStartupMessages({
  library(broom)
  library(data.table)
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
assert_saf_files(saf_paths)

pubsaf_dir <- saf_paths$pubsaf_dir
supp_dir <- saf_paths$supp_dir
pollution_dir <- file.path("output", "zip_pollution")
community_dir <- file.path("data", "processed", "community")
out_dir <- file.path("output", "diagnostics", "timevarying_pollution_models")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

analysis_start_year <- 2006L
analysis_end_year <- 2023L
target_organs <- c("HR", "KI", "LI", "LU")
max_intervals_per_organ <- as.integer(Sys.getenv("TV_DIAGNOSTIC_MAX_INTERVALS", "500000"))
include_center <- Sys.getenv("TV_DIAGNOSTIC_CENTER_FE", "FALSE") == "TRUE"

log_msg <- function(...) {
  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
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
  if (nzchar(detail)) message("\n", attr(pb, "label"), " ", value, "/", attr(pb, "total"), " | ", detail)
  invisible(NULL)
}

close_progress <- function(pb) {
  if (!is.null(pb)) close(pb)
}

clean_zip <- function(x) {
  z <- str_extract(as.character(x), "[0-9]+")
  z <- if_else(nchar(z) > 5L, substr(z, 1, 5), z)
  if_else(is.na(z), NA_character_, str_pad(z, 5, side = "left", pad = "0"))
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
    mutate(pm25_5ug = pm25_ug_m3 / 5, o3_10ppb = o3_ppb / 10, no2_10unit = no2 / 10)
}

candidate_files <- tribble(
  ~path,
  file.path(pubsaf_dir, "cand_kipa.sas7bdat"),
  file.path(pubsaf_dir, "cand_liin.sas7bdat"),
  file.path(pubsaf_dir, "cand_thor.sas7bdat")
)

candidate_cols <- c(
  "PERS_ID", "PX_ID", "WL_ORG", "CAN_LISTING_DT", "CAN_ACTIVATE_DT", "CAN_SOURCE",
  "CAN_REM_DT", "CAN_REM_CD", "CAN_DEATH_DT", "CAN_ENDWLFU", "CAN_AGE_AT_LISTING",
  "CAN_GENDER", "CAN_RACE_SRTR", "CAN_ETHNICITY_SRTR", "CAN_BMI", "CAN_LISTING_CTR_CD"
)

log_msg("Reading inputs")
candidate_zip <- read_sas(file.path(supp_dir, "canzip2506.sas7bdat")) %>%
  transmute(PERS_ID, PX_ID, candidate_zip = clean_zip(CAN_PERM_ZIP)) %>%
  distinct(PERS_ID, PX_ID, .keep_all = TRUE)

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

candidate_data <- candidate_files %>%
  rowwise() %>%
  mutate(data = list(read_sas(path, col_select = any_of(candidate_cols)))) %>%
  ungroup() %>%
  select(data) %>%
  tidyr::unnest(data)

baseline_exposure <- pollution_annual %>%
  transmute(
    candidate_zip = zip,
    prior_exposure_year = year,
    pm25_prior_1y_5ug = pm25_5ug,
    o3_prior_1y_10ppb = o3_10ppb,
    no2_prior_1y_10unit = no2_10unit
  )

log_msg("Constructing cohort")
cohort <- candidate_data %>%
  filter(WL_ORG %in% target_organs) %>%
  left_join(candidate_zip, by = c("PERS_ID", "PX_ID")) %>%
  mutate(
    waitlist_row_id = row_number(),
    index_date = coalesce(CAN_ACTIVATE_DT, CAN_LISTING_DT),
    index_year = as.integer(format(index_date, "%Y")),
    prior_exposure_year = index_year - 1L,
    adverse_waitlist_event = as.integer(CAN_REM_CD %in% c(8, 13)),
    transplant_event = as.integer(CAN_REM_CD %in% c(4, 14, 15, 18, 19, 21, 22, 23)),
    end_date = coalesce(CAN_REM_DT, CAN_DEATH_DT, CAN_ENDWLFU),
    followup_days = as.numeric(end_date - index_date),
    age = CAN_AGE_AT_LISTING,
    sex = as_factor_missing(if_else(CAN_GENDER %in% c("M", "F"), CAN_GENDER, NA_character_)),
    race_srtr = as_factor_missing(CAN_RACE_SRTR),
    ethnicity_srtr = as_factor_missing(CAN_ETHNICITY_SRTR),
    bmi = CAN_BMI,
    listing_center = as_factor_missing(CAN_LISTING_CTR_CD),
    index_year_factor = factor(index_year)
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
  mutate(followup_days = pmax(followup_days, 0.5)) %>%
  left_join(baseline_exposure, by = c("candidate_zip", "prior_exposure_year")) %>%
  left_join(community, by = c("candidate_zip" = "zip", "index_year" = "analysis_year")) %>%
  filter(!is.na(pm25_prior_1y_5ug), !is.na(o3_prior_1y_10ppb), !is.na(no2_prior_1y_10unit))

log_msg("Constructing annual waitlist intervals")
intervals <- as.data.table(cohort)
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
    group_by(WL_ORG) %>%
    summarise(
      n_intervals = n(),
      n_registrations = n_distinct(waitlist_row_id),
      adverse_events = sum(adverse_event_interval),
      sampled_for_models = n_intervals > max_intervals_per_organ,
      .groups = "drop"
    ),
  file.path(out_dir, "timevarying_diagnostic_interval_summary.csv")
)

pollutants <- tribble(
  ~pollutant, ~term,
  "PM2.5", "pm25_waitlist_cum_5ug",
  "Ozone", "o3_waitlist_cum_10ppb",
  "NO2", "no2_waitlist_cum_10unit"
)

base_terms <- c(
  "ns(age, df = 4)", "sex", "race_srtr", "ethnicity_srtr", "ns(bmi, df = 3)",
  "ns(median_household_income, df = 3)", "ns(pct_poverty, df = 3)",
  "ns(pct_bachelor_plus, df = 3)", "ns(pct_unemployed, df = 3)",
  "ns(pct_no_vehicle, df = 3)", "ns(pct_black, df = 3)", "ns(pct_hispanic, df = 3)",
  "strata(index_year_factor)", "cluster(PERS_ID)"
)
if (include_center) base_terms <- c(base_terms, "strata(listing_center)")

extract_formula_vars <- function(terms) {
  vars <- unique(unlist(str_extract_all(terms, "[A-Za-z][A-Za-z0-9_]*")))
  setdiff(vars, c("ns", "df", "strata", "cluster"))
}

wald_test_terms <- function(fit, exposure_term) {
  idx <- grep(paste0("(^|ns\\()", exposure_term), names(coef(fit)))
  if (length(idx) == 0L) return(tibble(wald_chisq = NA_real_, df = NA_integer_, p_value = NA_real_))
  beta <- coef(fit)[idx]
  v <- vcov(fit)[idx, idx, drop = FALSE]
  stat <- as.numeric(t(beta) %*% solve(v, beta))
  tibble(wald_chisq = stat, df = length(idx), p_value = pchisq(stat, df = length(idx), lower.tail = FALSE))
}

fit_spec <- function(org, model_type, exposure_form, terms_tbl) {
  dat <- as_tibble(intervals[WL_ORG == org])
  if (nrow(dat) > max_intervals_per_organ) {
    set.seed(20260528 + match(org, target_organs))
    events <- dat %>% filter(adverse_event_interval == 1L)
    nonevents <- dat %>% filter(adverse_event_interval == 0L)
    n_nonevents <- max(max_intervals_per_organ - nrow(events), 0L)
    dat <- bind_rows(events, slice_sample(nonevents, n = min(n_nonevents, nrow(nonevents)))) %>%
      arrange(waitlist_row_id, tstart)
  }

  exposure_terms <- if (exposure_form == "linear") terms_tbl$term else paste0("ns(", terms_tbl$term, ", df = 3)")
  rhs <- c(exposure_terms, base_terms)
  vars_needed <- intersect(extract_formula_vars(rhs), names(dat))
  model_dat <- dat %>%
    filter(complete.cases(across(all_of(vars_needed)))) %>%
    droplevels()
  if (nrow(model_dat) == 0L || sum(model_dat$adverse_event_interval) < 100L) return(NULL)

  form <- as.formula(paste("Surv(tstart, tstop, adverse_event_interval) ~", paste(rhs, collapse = " + ")))
  model_name <- paste(model_type, exposure_form, sep = "_")
  log_msg("Fitting ", org, " ", model_name, " n=", nrow(model_dat), " adverse=", sum(model_dat$adverse_event_interval))
  fit <- coxph(form, data = model_dat, ties = "efron", robust = TRUE, x = FALSE, y = FALSE)

  linear <- tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(term %in% terms_tbl$term) %>%
    mutate(
      organ = org,
      model_type = model_type,
      exposure_form = exposure_form,
      pollutant = recode(term, !!!setNames(terms_tbl$pollutant, terms_tbl$term)),
      n_intervals = nrow(model_dat),
      n_registrations = n_distinct(model_dat$waitlist_row_id),
      adverse_events = sum(model_dat$adverse_event_interval),
      sampled = nrow(dat) < nrow(as_tibble(intervals[WL_ORG == org])),
      concordance = unname(summary(fit)$concordance[1]),
      .before = 1
    )

  contrasts <- bind_rows(lapply(seq_len(nrow(terms_tbl)), function(i) {
    exposure_term <- terms_tbl$term[[i]]
    q <- quantile(model_dat[[exposure_term]], probs = c(0.25, 0.75), na.rm = TRUE)
    base_row <- model_dat[which.min(abs(model_dat[[exposure_term]] - median(model_dat[[exposure_term]], na.rm = TRUE))), , drop = FALSE]
    low_row <- base_row
    high_row <- base_row
    low_row[[exposure_term]] <- q[[1]]
    high_row[[exposure_term]] <- q[[2]]
    pred <- predict(fit, newdata = bind_rows(low_row, high_row), type = "lp", se.fit = TRUE)
    log_hr <- pred$fit[[2]] - pred$fit[[1]]
    se <- sqrt(pred$se.fit[[1]]^2 + pred$se.fit[[2]]^2)
    tibble(
      organ = org,
      model_type = model_type,
      exposure_form = exposure_form,
      pollutant = terms_tbl$pollutant[[i]],
      contrast = "p75_vs_p25",
      p25 = as.numeric(q[[1]]),
      p75 = as.numeric(q[[2]]),
      hazard_ratio = exp(log_hr),
      conf_low = exp(log_hr - 1.96 * se),
      conf_high = exp(log_hr + 1.96 * se),
      n_intervals = nrow(model_dat),
      n_registrations = n_distinct(model_dat$waitlist_row_id),
      adverse_events = sum(model_dat$adverse_event_interval),
      sampled = nrow(dat) < nrow(as_tibble(intervals[WL_ORG == org])),
      concordance = unname(summary(fit)$concordance[1])
    ) %>%
      bind_cols(wald_test_terms(fit, exposure_term))
  }))

  list(linear = linear, contrasts = contrasts)
}

all_results <- list()
total_specs <- length(target_organs) * (nrow(pollutants) * 2L + 2L)
pb <- make_progress(total_specs, "Time-varying diagnostic model grid")
progress_i <- 0L
for (org in target_organs) {
  for (i in seq_len(nrow(pollutants))) {
    terms_tbl <- pollutants[i, ]
    progress_i <- progress_i + 1L
    tick_progress(pb, progress_i, paste(org, terms_tbl$pollutant, "single linear"))
    all_results[[length(all_results) + 1L]] <- fit_spec(org, paste0("single_", str_to_lower(str_replace_all(terms_tbl$pollutant, "[^A-Za-z0-9]+", ""))), "linear", terms_tbl)
    progress_i <- progress_i + 1L
    tick_progress(pb, progress_i, paste(org, terms_tbl$pollutant, "single spline"))
    all_results[[length(all_results) + 1L]] <- fit_spec(org, paste0("single_", str_to_lower(str_replace_all(terms_tbl$pollutant, "[^A-Za-z0-9]+", ""))), "spline", terms_tbl)
  }
  progress_i <- progress_i + 1L
  tick_progress(pb, progress_i, paste(org, "combined linear"))
  all_results[[length(all_results) + 1L]] <- fit_spec(org, "combined", "linear", pollutants)
  progress_i <- progress_i + 1L
  tick_progress(pb, progress_i, paste(org, "combined spline"))
  all_results[[length(all_results) + 1L]] <- fit_spec(org, "combined", "spline", pollutants)
}
close_progress(pb)

write_csv(bind_rows(lapply(all_results, `[[`, "linear")), file.path(out_dir, "timevarying_single_vs_combined_linear_terms.csv"))
write_csv(bind_rows(lapply(all_results, `[[`, "contrasts")), file.path(out_dir, "timevarying_single_vs_combined_spline_iqr_contrasts.csv"))

log_msg("Wrote time-varying pollution diagnostics to ", out_dir)
