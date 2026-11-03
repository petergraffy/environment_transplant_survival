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
  library(readr)
  library(survival)
  library(tibble)
})

in_path <- file.path(
  "output",
  "primary_waitlist_period_pollution_cox",
  "primary_waitlist_period_pollution_analysis_dataset.csv.gz"
)
out_dir <- file.path("output", "primary_waitlist_period_pollution_subgroup_cox")
model_dir <- file.path(out_dir, "model_results")
dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)

organ_labels <- c(HR = "Heart", KI = "Kidney", LI = "Liver", LU = "Lung")
target_organs <- c("HR", "KI", "LI", "LU")

exposure_specs <- list(
  pm25 = list(
    term = "pm25_waitlist_5ug",
    raw_col = "pm25_waitlist_ug_m3",
    label = "PM2.5 waitlist-period mean per 5 ug/m3",
    followup = "followup_days_pm25",
    event = "event_pm25",
    end_year = 2023L
  ),
  o3 = list(
    term = "o3_waitlist_10ppb",
    raw_col = "o3_waitlist_ppb",
    label = "Ozone waitlist-period mean per 10 ppb",
    followup = "followup_days_o3",
    event = "event_o3",
    end_year = 2023L
  ),
  no2 = list(
    term = "no2_waitlist_10unit",
    raw_col = "no2_waitlist",
    label = "NO2 waitlist-period mean per 10 ppb",
    followup = "followup_days_no2",
    event = "event_no2",
    end_year = 2025L
  )
)

subgroup_specs <- list(
  age_group = list(
    label = "Age group",
    adjust = c("sex", "race_group", "listing_year", "organ_score", "strata(listing_center)")
  ),
  sex = list(
    label = "Sex",
    adjust = c("age", "race_group", "listing_year", "organ_score", "strata(listing_center)")
  ),
  race_group = list(
    label = "Race",
    adjust = c("age", "sex", "listing_year", "organ_score", "strata(listing_center)")
  )
)

organ_adjust_terms <- function(org, subgroup_spec) {
  terms <- subgroup_spec$adjust
  if (org == "KI") {
    terms <- c(
      terms[terms != "organ_score"],
      "kidney_no_dialysis_time",
      "kidney_dialysis_years",
      "kidney_diabetes"
    )
  }
  unique(terms)
}

adjust_vars_from_terms <- function(terms) {
  vars <- terms[!grepl("^strata\\(", terms)]
  gsub("^strata\\((.*)\\)$", "\\1", vars)
}

log_msg <- function(...) {
  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
  flush.console()
}

safe_error <- function(x) {
  if (is.null(x) || length(x) == 0L) return(NA_character_)
  as.character(x[[1]])
}

format_race_group <- function(x) {
  case_when(
    x == "WHITE" ~ "White",
    x == "BLACK" ~ "Black",
    x == "ASIAN" ~ "Asian",
    TRUE ~ "Other/Unknown"
  )
}

fit_one_subgroup <- function(dat, org, exposure_name, exposure_spec, subgroup_var, subgroup_level, subgroup_spec) {
  adjustment_terms <- organ_adjust_terms(org, subgroup_spec)
  adjustment_vars <- adjust_vars_from_terms(adjustment_terms)
  result_path <- file.path(
    model_dir,
    paste0(tolower(org), "_", exposure_name, "_", subgroup_var, "_", make.names(subgroup_level), ".csv")
  )
  vars_needed <- c(
    exposure_spec$followup,
    exposure_spec$event,
    exposure_spec$term,
    subgroup_var,
    "listing_center",
    "PERS_ID",
    "listing_year",
    adjustment_vars
  )
  if ("age" %in% subgroup_spec$adjust) vars_needed <- c(vars_needed, "age")
  if ("sex" %in% subgroup_spec$adjust) vars_needed <- c(vars_needed, "sex")
  if ("race_group" %in% subgroup_spec$adjust) vars_needed <- c(vars_needed, "race_group")

  model_dat <- dat %>%
    filter(
      WL_ORG == org,
      listing_year_int <= exposure_spec$end_year,
      .data[[subgroup_var]] == subgroup_level,
      complete.cases(across(all_of(unique(vars_needed)))),
      .data[[exposure_spec$followup]] > 0
    ) %>%
    mutate(
      .followup_days = pmax(.data[[exposure_spec$followup]], 0.5),
      .event = .data[[exposure_spec$event]]
    ) %>%
    droplevels()

  source_n <- nrow(model_dat)
  source_events <- sum(model_dat$.event, na.rm = TRUE)
  source_centers <- n_distinct(model_dat$listing_center)
  if (source_n == 0L || source_events < 50L || source_centers < 2L) {
    result <- tibble(
      organ = org,
      organ_label = recode(org, !!!organ_labels),
      exposure = exposure_name,
      exposure_label = exposure_spec$label,
      subgroup = subgroup_var,
      subgroup_label = subgroup_spec$label,
      subgroup_level = subgroup_level,
      exposure_data_end_year = exposure_spec$end_year,
      n = source_n,
      people = n_distinct(model_dat$PERS_ID),
      centers = source_centers,
      adverse_events = source_events,
      hazard_ratio = NA_real_,
      conf_low = NA_real_,
      conf_high = NA_real_,
      p_value = NA_real_,
      concordance = NA_real_,
      adjustment_set = paste(adjustment_terms, collapse = " + "),
      variance_estimator = "model_based",
      center_adjustment = "stratified_baseline_hazard_by_listing_center",
      error = "Insufficient sample/events/centers"
    )
    write_csv(result, result_path)
    return(result)
  }

  form <- as.formula(paste(
    "Surv(.followup_days, .event) ~",
    paste(c(exposure_spec$term, adjustment_terms), collapse = " + ")
  ))

  fit_error <- NULL
  fit <- tryCatch(
    coxph(form, data = model_dat, ties = "efron", x = FALSE, y = FALSE),
    error = function(e) {
      fit_error <<- conditionMessage(e)
      NULL
    }
  )

  if (is.null(fit)) {
    result <- tibble(
      organ = org,
      organ_label = recode(org, !!!organ_labels),
      exposure = exposure_name,
      exposure_label = exposure_spec$label,
      subgroup = subgroup_var,
      subgroup_label = subgroup_spec$label,
      subgroup_level = subgroup_level,
      exposure_data_end_year = exposure_spec$end_year,
      n = source_n,
      people = n_distinct(model_dat$PERS_ID),
      centers = source_centers,
      adverse_events = source_events,
      hazard_ratio = NA_real_,
      conf_low = NA_real_,
      conf_high = NA_real_,
      p_value = NA_real_,
      concordance = NA_real_,
      adjustment_set = paste(adjustment_terms, collapse = " + "),
      variance_estimator = "model_based",
      center_adjustment = "stratified_baseline_hazard_by_listing_center",
      error = safe_error(fit_error)
    )
  } else {
    result <- tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
      filter(term == exposure_spec$term) %>%
      transmute(
        organ = org,
        organ_label = recode(org, !!!organ_labels),
        exposure = exposure_name,
        exposure_label = exposure_spec$label,
        subgroup = subgroup_var,
        subgroup_label = subgroup_spec$label,
        subgroup_level = subgroup_level,
        exposure_data_end_year = exposure_spec$end_year,
        n = source_n,
        people = n_distinct(model_dat$PERS_ID),
        centers = source_centers,
        adverse_events = source_events,
        hazard_ratio = estimate,
        conf_low = conf.low,
        conf_high = conf.high,
        p_value = p.value,
        concordance = unname(summary(fit)$concordance[1]),
        adjustment_set = paste(adjustment_terms, collapse = " + "),
        variance_estimator = "model_based",
        center_adjustment = "stratified_baseline_hazard_by_listing_center",
        error = NA_character_
      )
  }

  write_csv(result, result_path)
  result
}

fit_interaction <- function(dat, org, exposure_name, exposure_spec, subgroup_var, subgroup_spec) {
  adjustment_terms <- organ_adjust_terms(org, subgroup_spec)
  adjustment_vars <- adjust_vars_from_terms(adjustment_terms)
  vars_needed <- c(
    exposure_spec$followup,
    exposure_spec$event,
    exposure_spec$term,
    subgroup_var,
    "listing_center",
    "PERS_ID",
    "listing_year",
    adjustment_vars
  )
  if ("age" %in% subgroup_spec$adjust) vars_needed <- c(vars_needed, "age")
  if ("sex" %in% subgroup_spec$adjust) vars_needed <- c(vars_needed, "sex")
  if ("race_group" %in% subgroup_spec$adjust) vars_needed <- c(vars_needed, "race_group")

  model_dat <- dat %>%
    filter(
      WL_ORG == org,
      listing_year_int <= exposure_spec$end_year,
      complete.cases(across(all_of(unique(vars_needed)))),
      .data[[exposure_spec$followup]] > 0
    ) %>%
    mutate(
      .followup_days = pmax(.data[[exposure_spec$followup]], 0.5),
      .event = .data[[exposure_spec$event]]
    ) %>%
    droplevels()

  source_n <- nrow(model_dat)
  source_events <- sum(model_dat$.event, na.rm = TRUE)
  source_centers <- n_distinct(model_dat$listing_center)
  if (source_n == 0L || source_events < 50L || source_centers < 2L || n_distinct(model_dat[[subgroup_var]]) < 2L) {
    return(tibble(
      organ = org,
      organ_label = recode(org, !!!organ_labels),
      exposure = exposure_name,
      subgroup = subgroup_var,
      subgroup_label = subgroup_spec$label,
      exposure_data_end_year = exposure_spec$end_year,
      n = source_n,
      people = n_distinct(model_dat$PERS_ID),
      centers = source_centers,
      adverse_events = source_events,
      subgroup_levels = n_distinct(model_dat[[subgroup_var]]),
      interaction_df = NA_real_,
      interaction_p_value = NA_real_,
      error = "Insufficient sample/events/centers/subgroup levels"
    ))
  }

  base_form <- as.formula(paste(
    "Surv(.followup_days, .event) ~",
    paste(c(exposure_spec$term, subgroup_var, adjustment_terms), collapse = " + ")
  ))
  int_form <- as.formula(paste(
    "Surv(.followup_days, .event) ~",
    paste(c(paste0(exposure_spec$term, " * ", subgroup_var), adjustment_terms), collapse = " + ")
  ))

  fit_error <- NULL
  base_fit <- tryCatch(
    coxph(base_form, data = model_dat, ties = "efron", x = FALSE, y = FALSE),
    error = function(e) {
      fit_error <<- conditionMessage(e)
      NULL
    }
  )
  int_fit <- if (!is.null(base_fit)) {
    tryCatch(
      coxph(int_form, data = model_dat, ties = "efron", x = FALSE, y = FALSE),
      error = function(e) {
        fit_error <<- conditionMessage(e)
        NULL
      }
    )
  } else {
    NULL
  }

  if (is.null(base_fit) || is.null(int_fit)) {
    return(tibble(
      organ = org,
      organ_label = recode(org, !!!organ_labels),
      exposure = exposure_name,
      subgroup = subgroup_var,
      subgroup_label = subgroup_spec$label,
      exposure_data_end_year = exposure_spec$end_year,
      n = source_n,
      people = n_distinct(model_dat$PERS_ID),
      centers = source_centers,
      adverse_events = source_events,
      subgroup_levels = n_distinct(model_dat[[subgroup_var]]),
      interaction_df = NA_real_,
      interaction_p_value = NA_real_,
      error = safe_error(fit_error)
    ))
  }

  lrt <- anova(base_fit, int_fit, test = "LRT")
  tibble(
    organ = org,
    organ_label = recode(org, !!!organ_labels),
    exposure = exposure_name,
    subgroup = subgroup_var,
    subgroup_label = subgroup_spec$label,
    exposure_data_end_year = exposure_spec$end_year,
    n = source_n,
    people = n_distinct(model_dat$PERS_ID),
    centers = source_centers,
    adverse_events = source_events,
    subgroup_levels = n_distinct(model_dat[[subgroup_var]]),
    interaction_df = lrt$Df[2],
    interaction_p_value = lrt$`Pr(>|Chi|)`[2],
    error = NA_character_
  )
}

log_msg("Reading primary waitlist-period pollution analysis dataset")
analysis_dat <- read_csv(in_path, show_col_types = FALSE) %>%
  mutate(
    pm25_waitlist_5ug = pm25_waitlist_ug_m3 / 5,
    o3_waitlist_10ppb = o3_waitlist_ppb / 10,
    no2_waitlist_10unit = no2_waitlist / 10,
    listing_year_int = as.integer(as.character(listing_year)),
    listing_year = factor(listing_year_int),
    sex = factor(sex),
    race_group = factor(format_race_group(race), levels = c("White", "Black", "Asian", "Other/Unknown")),
    age_group = cut(
      age,
      breaks = c(-Inf, 17, 39, 59, Inf),
      labels = c("<18", "18-39", "40-59", "60+"),
      right = TRUE
    ),
    age_group = factor(age_group, levels = c("<18", "18-39", "40-59", "60+")),
    listing_center = factor(listing_center)
  )

subgroup_results <- list()
interaction_results <- list()
counter <- 0L
total <- length(target_organs) * length(exposure_specs) * length(subgroup_specs)

for (org in target_organs) {
  for (exposure_name in names(exposure_specs)) {
    exposure_spec <- exposure_specs[[exposure_name]]
    for (subgroup_var in names(subgroup_specs)) {
      subgroup_spec <- subgroup_specs[[subgroup_var]]
      counter <- counter + 1L
      log_msg("Subgroup Cox ", counter, "/", total, " | ", org, " | ", exposure_name, " | ", subgroup_var)
      subgroup_levels <- levels(droplevels(analysis_dat[[subgroup_var]]))
      subgroup_results[[paste(org, exposure_name, subgroup_var, sep = "_")]] <- bind_rows(lapply(subgroup_levels, function(level) {
        fit_one_subgroup(analysis_dat, org, exposure_name, exposure_spec, subgroup_var, level, subgroup_spec)
      }))
      interaction_results[[paste(org, exposure_name, subgroup_var, sep = "_")]] <- fit_interaction(
        analysis_dat,
        org,
        exposure_name,
        exposure_spec,
        subgroup_var,
        subgroup_spec
      )
      invisible(gc())
    }
  }
}

subgroup_results <- bind_rows(subgroup_results)
interaction_results <- bind_rows(interaction_results)

write_csv(subgroup_results, file.path(out_dir, "primary_waitlist_period_pollution_subgroup_cox_results.csv"))
write_csv(interaction_results, file.path(out_dir, "primary_waitlist_period_pollution_subgroup_interactions.csv"))

subgroup_summary <- subgroup_results %>%
  group_by(organ, organ_label, exposure, subgroup, subgroup_label) %>%
  summarise(
    modeled_levels = sum(!is.na(hazard_ratio)),
    skipped_levels = sum(is.na(hazard_ratio)),
    min_events = min(adverse_events, na.rm = TRUE),
    max_events = max(adverse_events, na.rm = TRUE),
    .groups = "drop"
  )
write_csv(subgroup_summary, file.path(out_dir, "primary_waitlist_period_pollution_subgroup_summary.csv"))

log_msg("Wrote subgroup Cox outputs to ", out_dir)
