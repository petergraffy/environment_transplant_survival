#!/usr/bin/env Rscript

runtime_source <- file.path("code", "r_runtime.R")
if (!file.exists(runtime_source)) runtime_source <- "r_runtime.R"
if (file.exists(runtime_source)) {
  source(runtime_source)
  ensure_user_library()
}

suppressPackageStartupMessages({
  library(cmprsk)
  library(data.table)
  library(dplyr)
  library(readr)
  library(tibble)
})

out_dir <- file.path("output", "formal_waitlist_environment_timevarying")
fg_dir <- file.path(out_dir, "reduced_finegray_heat_humidity_results")
dir.create(fg_dir, recursive = TRUE, showWarnings = FALSE)

target_organs <- c("HR", "KI", "LI", "LU")
organ_labels <- c(HR = "Heart", KI = "Kidney", LI = "Liver", LU = "Lung")
finegray_max_n <- as.integer(Sys.getenv("FINEGRAY_HEAT_MAX_N", "50000"))

log_msg <- function(...) {
  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
  flush.console()
}

log_msg("Reading static formal dataset and month/day-weighted heat exposures")
base <- read_csv(
  file.path("output", "formal_waitlist_environment", "formal_waitlist_environment_analysis_dataset.csv.gz"),
  show_col_types = FALSE
) %>%
  select(
    waitlist_row_id, PERS_ID, WL_ORG, followup_days, event_type,
    age, sex, race, index_year_centered, organ_score
  )

heat <- read_csv(
  file.path("output", "waitlist_environment_exposure", "waitlist_environment_exposures.csv.gz"),
  show_col_types = FALSE
) %>%
  select(waitlist_row_id, tmax_c_waitlist_mean, rmax_pct_waitlist_mean)

dat <- base %>%
  inner_join(heat, by = "waitlist_row_id") %>%
  mutate(
    sex = factor(sex),
    race = factor(race),
    tmax_month_weighted_5c = tmax_c_waitlist_mean / 5,
    rmax_month_weighted_10pct = rmax_pct_waitlist_mean / 10
  ) %>%
  as.data.table()

fit_reduced_fg <- function(model_dat, org) {
  model_dat <- model_dat[
    WL_ORG == org &
      complete.cases(
        followup_days, event_type, age, sex, race, index_year_centered,
        organ_score, tmax_month_weighted_5c, rmax_month_weighted_10pct
      )
  ]
  if (nrow(model_dat) == 0L || sum(model_dat$event_type == "adverse") < 100L) return(tibble())

  sampled <- FALSE
  source_n <- nrow(model_dat)
  source_adverse <- sum(model_dat$event_type == "adverse")
  if (source_n > finegray_max_n) {
    set.seed(20260603 + match(org, target_organs))
    sampled <- TRUE
    model_dat <- as_tibble(model_dat) %>%
      group_by(event_type) %>%
      slice_sample(prop = min(1, finegray_max_n / source_n)) %>%
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
    model = "reduced_finegray_heat_humidity_month_day_weighted",
    endpoint = "death_or_deterioration_delist_subdistribution",
    term = keep,
    exposure = recode(
      keep,
      tmax_month_weighted_5c = "Month/day-weighted Tmax per 5 C",
      rmax_month_weighted_10pct = "Month/day-weighted maximum relative humidity per 10 pct"
    ),
    source_n = source_n,
    source_adverse_events = source_adverse,
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

paths <- character()
for (org in target_organs) {
  path <- file.path(fg_dir, paste0(tolower(org), "_heat_humidity.csv"))
  paths <- c(paths, path)
  if (file.exists(path)) {
    log_msg("Skipping existing reduced Fine-Gray result ", org)
    next
  }
  result <- fit_reduced_fg(dat, org)
  write_csv(result, path)
  rm(result)
  invisible(gc())
}

results <- bind_rows(lapply(paths[file.exists(paths)], read_csv, show_col_types = FALSE))
write_csv(results, file.path(out_dir, "reduced_finegray_heat_humidity_results.csv"))
log_msg("Wrote reduced Fine-Gray heat/humidity sensitivity to ", out_dir)
