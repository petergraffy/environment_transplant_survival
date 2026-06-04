#!/usr/bin/env Rscript

runtime_source <- file.path("code", "r_runtime.R")
if (!file.exists(runtime_source)) runtime_source <- "r_runtime.R"
if (file.exists(runtime_source)) {
  source(runtime_source)
  ensure_user_library()
}

suppressPackageStartupMessages({
  library(broom)
  library(cmprsk)
  library(data.table)
  library(dplyr)
  library(readr)
  library(survival)
  library(tibble)
})

out_dir <- file.path("output", "formal_waitlist_environment_timevarying")
model_result_dir <- file.path(out_dir, "timevarying_cox_model_results")
dir.create(model_result_dir, recursive = TRUE, showWarnings = FALSE)

target_organs <- c("HR", "KI", "LI", "LU")
organ_labels <- c(HR = "Heart", KI = "Kidney", LI = "Liver", LU = "Lung")

log_msg <- function(...) {
  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
  flush.console()
}

log_msg("Reading saved time-varying interval dataset")
tv <- read_csv(file.path(out_dir, "timevarying_interval_analysis_dataset.csv.gz"), show_col_types = FALSE) %>%
  mutate(
    sex = factor(sex),
    race = factor(race),
    listing_center = factor(listing_center),
    tmax_interval_5c = tmax_interval_c / 5,
    rmax_interval_10pct = rmax_interval_pct / 10,
    pm25_interval_5ug = pm25_interval_ug_m3 / 5,
    o3_interval_10ppb = o3_interval_ppb / 10,
    no2_interval_10unit = no2_interval / 10
  ) %>%
  as.data.table()

exposure_specs <- list(
  heat_humidity = c("tmax_interval_5c", "rmax_interval_10pct"),
  pm25 = "pm25_interval_5ug",
  o3 = "o3_interval_10ppb",
  no2 = "no2_interval_10unit"
)

term_labels <- c(
  tmax_interval_5c = "Interval Tmax per 5 C",
  rmax_interval_10pct = "Interval maximum relative humidity per 10 pct",
  pm25_interval_5ug = "Interval PM2.5 per 5 ug/m3",
  o3_interval_10ppb = "Interval ozone per 10 ppb",
  no2_interval_10unit = "Interval NO2 per 10 units"
)

fit_tv_cox <- function(dat, org, model_name, exposure_terms) {
  vars_needed <- c("tstart", "tstop", "tv_adverse_event", "age_interval", "sex", "race", "index_year_centered", "organ_score_tv", "listing_center", exposure_terms)
  model_dat <- dat[complete.cases(dat[, ..vars_needed]) & tstop > tstart]
  if (nrow(model_dat) == 0L || sum(model_dat$tv_adverse_event) < 100L) return(tibble())
  log_msg("Time-varying Cox ", org, " | ", model_name, " intervals=", nrow(model_dat), " adverse=", sum(model_dat$tv_adverse_event))
  rhs <- c(exposure_terms, "age_interval", "sex", "race", "index_year_centered", "organ_score_tv", "strata(listing_center)", "cluster(PERS_ID)")
  form <- as.formula(paste("Surv(tstart, tstop, tv_adverse_event) ~", paste(rhs, collapse = " + ")))
  fit <- coxph(form, data = model_dat, ties = "efron", robust = TRUE, x = FALSE, y = FALSE)
  tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(term %in% exposure_terms) %>%
    transmute(
      organ = org,
      organ_label = recode(org, !!!organ_labels),
      model = model_name,
      endpoint = "time_varying_death_or_deterioration_delist_cause_specific",
      term,
      exposure = recode(term, !!!term_labels),
      intervals = nrow(model_dat),
      waitlist_rows = uniqueN(model_dat$waitlist_row_id),
      people = uniqueN(model_dat$PERS_ID),
      adverse_events = sum(model_dat$tv_adverse_event),
      hazard_ratio = estimate,
      conf_low = conf.low,
      conf_high = conf.high,
      p_value = p.value,
      concordance = unname(summary(fit)$concordance[1]),
      variance_estimator = "robust_patient_clustered",
      center_adjustment = "stratified_listing_center_baseline_hazard"
    )
}

tv_result_paths <- character()
for (org in target_organs) {
  dat <- tv[WL_ORG == org]
  for (model_name in names(exposure_specs)) {
    result_path <- file.path(model_result_dir, paste0(tolower(org), "_", model_name, ".csv"))
    tv_result_paths <- c(tv_result_paths, result_path)
    if (file.exists(result_path)) {
      log_msg("Skipping existing time-varying Cox result ", org, " | ", model_name)
      next
    }
    result <- fit_tv_cox(dat, org, model_name, exposure_specs[[model_name]])
    write_csv(result, result_path)
    rm(result)
    invisible(gc())
  }
}

tv_results <- bind_rows(lapply(tv_result_paths[file.exists(tv_result_paths)], read_csv, show_col_types = FALSE))
write_csv(tv_results, file.path(out_dir, "timevarying_cox_exposure_results.csv"))

log_msg("Fitting reduced Fine-Gray heat/humidity sensitivity")
fg_base <- tv[
  ,
  .(
    followup_days = max(tstop, na.rm = TRUE),
    event_observed = max(tv_adverse_event, na.rm = TRUE),
    age = first(age_interval),
    sex = first(sex),
    race = first(race),
    index_year_centered = first(index_year_centered),
    organ_score = first(organ_score_tv),
    tmax_month_weighted_5c = weighted.mean(tmax_interval_5c, pmax(tstop - tstart, 0.5), na.rm = TRUE),
    rmax_month_weighted_10pct = weighted.mean(rmax_interval_10pct, pmax(tstop - tstart, 0.5), na.rm = TRUE)
  ),
  by = .(waitlist_row_id, PERS_ID, WL_ORG)
]

static_events <- read_csv(
  file.path("output", "formal_waitlist_environment", "formal_waitlist_environment_analysis_dataset.csv.gz"),
  show_col_types = FALSE
) %>%
  select(waitlist_row_id, event_type) %>%
  as.data.table()
fg_base <- static_events[fg_base, on = "waitlist_row_id"]

fit_reduced_fg <- function(dat, org) {
  model_dat <- dat[
    WL_ORG == org &
      complete.cases(followup_days, event_type, age, sex, race, index_year_centered, organ_score, tmax_month_weighted_5c, rmax_month_weighted_10pct)
  ]
  if (nrow(model_dat) == 0L || sum(model_dat$event_type == "adverse") < 100L) return(tibble())
  finegray_max_n <- as.integer(Sys.getenv("FINEGRAY_HEAT_MAX_N", "100000"))
  sampled <- FALSE
  if (nrow(model_dat) > finegray_max_n) {
    set.seed(20260603 + match(org, target_organs))
    sampled <- TRUE
    model_dat <- as_tibble(model_dat) %>%
      group_by(event_type) %>%
      slice_sample(prop = min(1, finegray_max_n / nrow(model_dat))) %>%
      ungroup() %>%
      as.data.table()
  }
  log_msg("Reduced Fine-Gray heat/humidity ", org, " n=", nrow(model_dat), " adverse=", sum(model_dat$event_type == "adverse"))
  fstatus <- fcase(
    model_dat$event_type == "adverse", 1L,
    model_dat$event_type %in% c("transplant_or_improvement", "other_exit"), 2L,
    default = 0L
  )
  mm <- model.matrix(
    ~ tmax_month_weighted_5c + rmax_month_weighted_10pct + age + sex + race + index_year_centered + organ_score,
    data = model_dat
  )
  mm <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]
  fit <- crr(model_dat$followup_days, fstatus, cov1 = mm, failcode = 1, cencode = 0)
  coef <- fit$coef
  var <- diag(fit$var)
  keep <- intersect(c("tmax_month_weighted_5c", "rmax_month_weighted_10pct"), names(coef))
  tibble(
    organ = org,
    organ_label = recode(org, !!!organ_labels),
    model = "reduced_finegray_heat_humidity_month_weighted",
    endpoint = "death_or_deterioration_delist_subdistribution",
    term = keep,
    exposure = recode(
      keep,
      tmax_month_weighted_5c = "Month/day-weighted Tmax per 5 C",
      rmax_month_weighted_10pct = "Month/day-weighted maximum relative humidity per 10 pct"
    ),
    n = nrow(model_dat),
    adverse_events = sum(model_dat$event_type == "adverse"),
    sampled = sampled,
    finegray_max_n = finegray_max_n,
    subdistribution_hazard_ratio = exp(coef[keep]),
    conf_low = exp(coef[keep] - 1.96 * sqrt(var[match(keep, names(coef))])),
    conf_high = exp(coef[keep] + 1.96 * sqrt(var[match(keep, names(coef))])),
    p_value = 2 * pnorm(abs(coef[keep] / sqrt(var[match(keep, names(coef))])), lower.tail = FALSE),
    center_adjustment = "not_included_reduced_sensitivity",
    variance_estimator = "cmprsk_crr_model_based"
  )
}

fg_results <- bind_rows(lapply(target_organs, function(org) fit_reduced_fg(fg_base, org)))
write_csv(fg_results, file.path(out_dir, "reduced_finegray_heat_humidity_results.csv"))

write_csv(
  as_tibble(tv[, .(
    intervals = .N,
    waitlist_rows = uniqueN(waitlist_row_id),
    people = uniqueN(PERS_ID),
    adverse_events = sum(tv_adverse_event),
    complete_heat_humidity_intervals = sum(complete.cases(tmax_interval_5c, rmax_interval_10pct)),
    complete_pm25_intervals = sum(!is.na(pm25_interval_5ug)),
    complete_o3_intervals = sum(!is.na(o3_interval_10ppb)),
    complete_no2_intervals = sum(!is.na(no2_interval_10unit))
  ), by = WL_ORG]),
  file.path(out_dir, "timevarying_interval_summary.csv")
)

log_msg("Wrote time-varying model outputs to ", out_dir)
