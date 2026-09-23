#!/usr/bin/env Rscript

runtime_source <- file.path("code", "r_runtime.R")
if (!file.exists(runtime_source)) runtime_source <- "r_runtime.R"
if (file.exists(runtime_source)) {
  source(runtime_source)
  ensure_user_library()
}

source(file.path("code", "daily_pollution_inputs.R"))

suppressPackageStartupMessages({
  library(arrow)
  library(data.table)
  library(dplyr)
  library(fixest)
  library(ggplot2)
  library(haven)
  library(patchwork)
  library(readr)
  library(scales)
  library(stringr)
  library(tibble)
  library(tidyr)
})

saf_paths_source <- file.path("code", "saf_paths.R")
if (!file.exists(saf_paths_source)) saf_paths_source <- "saf_paths.R"
source(saf_paths_source)
saf_paths <- get_saf_paths(release = "q1_2026")
pubsaf_dir <- saf_paths$pubsaf_dir

analysis_path <- file.path(
  "output",
  "primary_waitlist_period_pollution_cox",
  "primary_waitlist_period_pollution_analysis_dataset.csv.gz"
)
community_path <- file.path("data", "processed", "community", "zcta_acs_community_covariates_2005_2023.csv.gz")
release_dir <- file.path("data", "release")
annual_pollution_dir <- file.path(release_dir, "air_pollution_zcta_parquet")
out_dir <- file.path("output", "timevarying_pollution_severity_models_with_listing_year")
fig_dir <- file.path("output", "figures", "pollution_severity_associations_with_listing_year")
cache_dir <- file.path("output", "timevarying_pollution_cox_svi", "cache")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

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

flag_yes <- function(x) {
  y <- str_to_upper(str_trim(as.character(x)))
  as.integer(y %in% c("1", "Y", "YES", "TRUE", "T"))
}

safe_log1p <- function(x) {
  log(pmax(as.numeric(x), 0) + 1)
}

median_impute_num <- function(x) {
  x <- as.numeric(x)
  med <- median(x, na.rm = TRUE)
  if (!is.finite(med)) med <- NA_real_
  fifelse(is.na(x), med, x)
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

parquet_files <- function(path) {
  list.files(path, pattern = "[.]parquet$", full.names = TRUE)
}

make_complete_acs_svi_proxy <- function(path) {
  community <- read_csv(path, show_col_types = FALSE) %>%
    mutate(zip = clean_zip(zip), analysis_year = as.integer(analysis_year))
  community_2023 <- community %>% filter(analysis_year == 2023L) %>% select(-analysis_year)
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
    mutate(across(all_of(vulnerability_vars), ~ if_else(is.na(.x), median(.x, na.rm = TRUE), .x))) %>%
    mutate(
      svi_poverty_rank = percent_rank(pct_poverty),
      svi_unemployed_rank = percent_rank(pct_unemployed),
      svi_no_vehicle_rank = percent_rank(pct_no_vehicle),
      svi_nonwhite_rank = percent_rank(pct_nonwhite),
      svi_low_income_rank = percent_rank(-median_household_income),
      svi_low_education_rank = percent_rank(-pct_bachelor_plus),
      zcta_svi_proxy = rowMeans(cbind(
        svi_poverty_rank,
        svi_unemployed_rank,
        svi_no_vehicle_rank,
        svi_nonwhite_rank,
        svi_low_income_rank,
        svi_low_education_rank
      ), na.rm = TRUE)
    ) %>%
    ungroup() %>%
    select(zip, analysis_year, zcta_svi_proxy)
}

daily_to_monthly_cache <- function(path, value_col, out_col, cache_file) {
  read_daily_pollution_aggregate(path, value_col, out_col, "monthly", cache_file)
}

read_exposure_tables <- function() {
  pm25_monthly <- daily_to_monthly_cache(
    file.path(release_dir, "lghap_pm25_zcta_daily_parquet"),
    "pm25_ug_m3",
    "pm25_interval_ug_m3",
    file.path(cache_dir, "pm25_daily_derived_monthly_zcta.csv.gz")
  )
  o3_monthly <- daily_to_monthly_cache(
    file.path(release_dir, "o3_zcta_daily_parquet"),
    "o3_ppb",
    "o3_interval_ppb",
    file.path(cache_dir, "o3_daily_derived_monthly_zcta.csv.gz")
  )
  no2_monthly <- open_dataset(parquet_files(file.path(release_dir, "no2_zcta_monthly_parquet"))) %>%
    transmute(zip = zip, year = year, month = month, no2_interval_ppb = no2_ppbv) %>%
    collect() %>%
    mutate(zip = clean_zip(zip), year = as.integer(year), month = as.integer(month))
  no2_annual <- read_parquet(file.path(annual_pollution_dir, "air_pollution_zcta_no2_annual_2005_2025.parquet")) %>%
    transmute(zip = clean_zip(zip), year = as.integer(year), no2_annual_ppb = no2)
  list(
    pm25 = as.data.table(pm25_monthly),
    o3 = as.data.table(o3_monthly),
    no2_monthly = as.data.table(no2_monthly),
    no2_annual = as.data.table(no2_annual)
  )
}

make_month_intervals <- function(cohort_dt, exposure_end_date) {
  intervals <- trim_monthly_cohort(cohort_dt)
  intervals[, analysis_end_date := pmin(observed_end_date, exposure_end_date)]
  intervals <- intervals[analysis_end_date >= index_date]
  out <- expand_monthly_cohort(intervals)
  out[, interval_start_date := pmax(index_date, month_start)]
  out[, interval_end_date := pmin(analysis_end_date, month_end)]
  out <- out[interval_end_date >= interval_start_date]
  out[, `:=`(
    year = as.integer(format(month_start, "%Y")),
    month = as.integer(format(month_start, "%m")),
    tstart = as.numeric(interval_start_date - index_date),
    tstop = as.numeric(interval_end_date - index_date) + 1
  )]
  out[, kidney_dialysis_years_tv := kidney_dialysis_years]
  out[WL_ORG == "KI" & kidney_no_dialysis_time == 0L & !is.na(kidney_dialysis_years), kidney_dialysis_years_tv := kidney_dialysis_years + tstart / 365.25]
  out[, organ_score_tv := organ_score]
  out[, interval_row_id := .I]
  out
}

make_time_updated_score_sources <- function() {
  log_msg("Reading time-updated liver MELD and heart US-CRS inputs")
  liver_updates <- read_sas(
    file.path(pubsaf_dir, "stathist_liin.sas7bdat"),
    col_select = any_of(c("PX_ID", "WL_ORG", "CANHX_BEGIN_DT", "CANHX_SRTR_LAB_MELD", "CANHX_OPTN_LAB_MELD"))
  ) %>%
    filter(WL_ORG == "LI", !is.na(PX_ID), !is.na(CANHX_BEGIN_DT)) %>%
    transmute(
      PX_ID,
      score_date = as.Date(CANHX_BEGIN_DT),
      liver_score_update = coalesce(CANHX_SRTR_LAB_MELD, CANHX_OPTN_LAB_MELD),
      liver_score_update = if_else(!is.na(liver_score_update) & liver_score_update > 1000, liver_score_update - 6200, liver_score_update)
    ) %>%
    filter(is.finite(liver_score_update)) %>%
    as.data.table()

  heart_updates <- read_sas(
    file.path(pubsaf_dir, "statjust_hr1a.sas7bdat"),
    col_select = any_of(c(
      "PX_ID", "WL_ORG", "CANHX_CHG_DT", "CANHX_RVAD_TYPE", "CANHX_ECMO",
      "CANHX_LVAD_TYPE", "CANHX_LAB_SERUM_CREAT", "CANHX_LAB_BILI",
      "CANHX_LAB_ALBUMIN", "CANHX_LAB_SODIUM", "CANHX_LAB_BNP"
    ))
  ) %>%
    filter(WL_ORG == "HR", !is.na(PX_ID), !is.na(CANHX_CHG_DT)) %>%
    transmute(
      PX_ID,
      score_date = as.Date(CANHX_CHG_DT),
      hr_short_mcs_update = as.integer(flag_yes(CANHX_ECMO) == 1L | flag_yes(CANHX_RVAD_TYPE) == 1L),
      hr_durable_lvad_update = flag_yes(CANHX_LVAD_TYPE),
      hr_creatinine_update = as.numeric(CANHX_LAB_SERUM_CREAT),
      hr_bilirubin_update = as.numeric(CANHX_LAB_BILI),
      hr_albumin_update = as.numeric(CANHX_LAB_ALBUMIN),
      hr_sodium_update = as.numeric(CANHX_LAB_SODIUM),
      hr_bnp_update = as.numeric(CANHX_LAB_BNP)
    ) %>%
    as.data.table()

  if (nrow(heart_updates) > 0L) {
    heart_updates[, `:=`(
      hr_creatinine_update = median_impute_num(hr_creatinine_update),
      hr_bilirubin_update = median_impute_num(hr_bilirubin_update),
      hr_albumin_update = median_impute_num(hr_albumin_update),
      hr_sodium_update = median_impute_num(hr_sodium_update),
      hr_bnp_update = median_impute_num(hr_bnp_update)
    )]
  }

  setkey(liver_updates, PX_ID, score_date)
  setkey(heart_updates, PX_ID, score_date)
  list(liver = liver_updates, heart = heart_updates)
}

add_time_updated_scores <- function(intervals, score_sources) {
  intervals[, age_interval := age + tstart / 365.25]

  li_intervals <- intervals[WL_ORG == "LI", .(interval_row_id, PX_ID, interval_start_date)]
  if (nrow(li_intervals) > 0L && nrow(score_sources$liver) > 0L) {
    setkey(li_intervals, PX_ID, interval_start_date)
    li_join <- score_sources$liver[li_intervals, on = c("PX_ID", "score_date" = "interval_start_date"), roll = Inf]
    intervals[li_join, organ_score_tv := i.liver_score_update, on = "interval_row_id"]
  }

  hr_intervals <- intervals[WL_ORG == "HR", .(interval_row_id, PX_ID, interval_start_date, age_interval, sex)]
  if (nrow(hr_intervals) > 0L && nrow(score_sources$heart) > 0L) {
    setkey(hr_intervals, PX_ID, interval_start_date)
    hr_join <- score_sources$heart[hr_intervals, on = c("PX_ID", "score_date" = "interval_start_date"), roll = Inf]
    hr_join[, hr_egfr_update := calc_egfr_2021(hr_creatinine_update, age_interval, sex)]
    hr_join[, heart_score_update := 1.02 * hr_short_mcs_update +
      0.55 * safe_log1p(hr_bilirubin_update) -
      0.01 * hr_egfr_update +
      0.40 * safe_log1p(hr_bnp_update) -
      0.63 * hr_albumin_update -
      0.07 * hr_sodium_update -
      1.12 * hr_durable_lvad_update]
    intervals[hr_join, organ_score_tv := i.heart_score_update, on = "interval_row_id"]
  }

  intervals[is.na(organ_score_tv), organ_score_tv := organ_score]
  intervals
}

severity_outcome <- function(org) {
  if (org == "KI") {
    return(list(term = "kidney_dialysis_years_tv", label = "Dialysis duration, years"))
  }
  if (org == "LI") return(list(term = "organ_score_tv", label = "MELD/PELD"))
  if (org == "HR") return(list(term = "organ_score_tv", label = "US-CRS proxy"))
  if (org == "LU") return(list(term = "organ_score_tv", label = "Lung urgency proxy"))
  stop("Unsupported organ: ", org, call. = FALSE)
}

fit_tv_severity_model <- function(intervals, org, pollutant, exposure_term, exposure_label) {
  outcome_spec <- severity_outcome(org)
  vars_needed <- c(outcome_spec$term, exposure_term, "age_interval", "sex", "race", "listing_year", "zcta_svi_proxy", "listing_center")
  model_dat <- intervals[
    WL_ORG == org &
      complete.cases(intervals[, ..vars_needed])
  ]
  if (nrow(model_dat) == 0L) return(tibble())

  log_msg("Time-varying severity model ", org, " | ", pollutant, " intervals=", nrow(model_dat), " episodes=", uniqueN(model_dat$waitlist_row_id))
  form <- as.formula(paste0(outcome_spec$term, " ~ ", exposure_term, " + age_interval + sex + race + listing_year + zcta_svi_proxy | listing_center"))
  fit <- feols(form, data = model_dat, notes = FALSE)
  ci <- as.numeric(confint(fit, parm = exposure_term))
  beta <- coef(fit)[[exposure_term]]
  coef_table <- coeftable(fit)
  p_col <- grep("^Pr\\(", colnames(coef_table), value = TRUE)[1]
  pval <- coef_table[exposure_term, p_col]

  tibble(
    model_type = "time_varying",
    organ = org,
    organ_label = recode(org, !!!organ_labels),
    outcome = outcome_spec$term,
    outcome_label = outcome_spec$label,
    outcome_family = "gaussian",
    pollutant = pollutant,
    exposure = exposure_label,
    intervals = nrow(model_dat),
    candidate_organ_episodes = uniqueN(model_dat$waitlist_row_id),
    people = uniqueN(model_dat$PERS_ID),
    centers = uniqueN(model_dat$listing_center),
    estimate = beta,
    conf_low = ci[[1]],
    conf_high = ci[[2]],
    p_value = pval,
    measure = "mean_difference",
    adjustment_set = "age + sex + race + listing_year + zcta_svi_proxy + listing_center fixed effects"
  )
}

plot_severity_associations <- function(baseline_results, tv_results) {
  continuous_baseline <- baseline_results %>%
    filter(measure == "mean_difference") %>%
    mutate(
      model_type = "baseline",
      intervals = NA_integer_,
      candidate_organ_episodes = n,
      people = NA_integer_,
      centers = NA_integer_
    ) %>%
    select(model_type, organ_label, outcome_label, pollutant, exposure, n = candidate_organ_episodes, estimate, conf_low, conf_high, p_value, measure)

  continuous_tv <- tv_results %>%
    transmute(model_type, organ_label, outcome_label, pollutant, exposure, n = candidate_organ_episodes, estimate, conf_low, conf_high, p_value, measure)

  plot_dat <- bind_rows(continuous_baseline, continuous_tv) %>%
    filter(pollutant %in% c("pm25", "o3", "no2")) %>%
    mutate(
      model_type = recode(model_type, baseline = "Baseline severity", time_varying = "Time-updated severity"),
      model_type = factor(model_type, levels = c("Baseline severity", "Time-updated severity")),
      pollutant = factor(pollutant, levels = c("pm25", "o3", "no2")),
      outcome_panel = case_when(
        organ_label == "Kidney" ~ "Kidney: dialysis duration",
        organ_label == "Liver" ~ "Liver: MELD/PELD",
        organ_label == "Heart" ~ "Heart: US-CRS proxy",
        organ_label == "Lung" ~ "Lung: urgency proxy",
        TRUE ~ outcome_label
      ),
      outcome_panel = factor(outcome_panel, levels = rev(c(
        "Kidney: dialysis duration",
        "Liver: MELD/PELD",
        "Heart: US-CRS proxy",
        "Lung: urgency proxy"
      )))
    )

  diabetes_dat <- baseline_results %>%
    filter(measure == "odds_ratio", outcome_label == "Diabetes at listing") %>%
    mutate(
      pollutant = factor(pollutant, levels = c("pm25", "o3", "no2")),
      outcome_panel = factor("Kidney: diabetes at listing")
    )

  pollutant_labels <- c(
    pm25 = 'PM[2.5]~"per 5 ug/m"^3',
    o3 = 'O[3]~"per 10 ppb"',
    no2 = 'NO[2]~"per 10 ppb"'
  )
  pollutant_colors <- c(pm25 = "#0072B2", o3 = "#009E73", no2 = "#D55E00")
  dodge <- position_dodge(width = 0.72)

  p_cont <- ggplot(plot_dat, aes(x = estimate, y = outcome_panel, color = pollutant, shape = model_type)) +
    geom_vline(xintercept = 0, color = "grey35", linewidth = 0.4) +
    geom_errorbar(aes(xmin = conf_low, xmax = conf_high), orientation = "y", width = 0.18, linewidth = 0.65, position = dodge) +
    geom_point(size = 3.2, position = dodge) +
    scale_color_manual(values = pollutant_colors, breaks = names(pollutant_labels), labels = parse(text = pollutant_labels)) +
    scale_shape_manual(values = c("Baseline severity" = 16, "Time-updated severity" = 17)) +
    labs(x = "Adjusted mean difference in severity measure", y = NULL, color = NULL, shape = NULL) +
    guides(
      color = guide_legend(nrow = 1, byrow = TRUE, order = 1),
      shape = guide_legend(nrow = 1, byrow = TRUE, order = 2)
    ) +
    theme_minimal(base_size = 17) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text = element_text(color = "grey15"),
      legend.position = "bottom",
      legend.box = "vertical",
      legend.justification = "center",
      legend.margin = margin(t = 6),
      plot.margin = margin(8, 12, 4, 8)
    )

  p_diab <- ggplot(diabetes_dat, aes(x = estimate, y = outcome_panel, color = pollutant)) +
    geom_vline(xintercept = 1, color = "grey35", linewidth = 0.4) +
    geom_errorbar(aes(xmin = conf_low, xmax = conf_high), orientation = "y", width = 0.18, linewidth = 0.65, position = dodge) +
    geom_point(size = 3.2, position = dodge) +
    scale_color_manual(values = pollutant_colors, breaks = names(pollutant_labels), labels = parse(text = pollutant_labels)) +
    scale_x_log10(labels = label_number(accuracy = 0.01)) +
    labs(x = "Adjusted odds ratio for diabetes at listing", y = NULL, color = NULL) +
    theme_minimal(base_size = 17) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text = element_text(color = "grey15"),
      legend.position = "none",
      plot.margin = margin(4, 12, 8, 8)
    )

  combined <- p_cont / p_diab + patchwork::plot_layout(heights = c(4, 1.2))

  png_path <- file.path(fig_dir, "pollution_severity_association_forest.png")
  pdf_path <- file.path(fig_dir, "pollution_severity_association_forest.pdf")
  ggsave(png_path, combined, width = 15.5, height = 9.2, dpi = 360, bg = "white")
  ggsave(pdf_path, combined, width = 15.5, height = 9.2, device = cairo_pdf, bg = "white")

  write_csv(
    tibble(
      figure = "pollution_severity_association_forest",
      path = normalizePath(c(png_path, pdf_path), winslash = "/", mustWork = FALSE)
    ),
    file.path(fig_dir, "pollution_severity_association_forest_manifest.csv")
  )
}

log_msg("Reading primary deduplicated cohort")
analysis_dat <- read_csv(analysis_path, show_col_types = FALSE) %>%
  mutate(
    candidate_zip = clean_zip(candidate_zip),
    index_date = as.Date(index_date),
    observed_end_date = as.Date(observed_end_date),
    listing_year_int = as.integer(as.character(listing_year)),
    listing_year = factor(listing_year_int),
    sex = factor(sex),
    race = factor(race),
    listing_center = factor(listing_center)
  ) %>%
  filter(WL_ORG %in% target_organs, observed_end_date >= index_date, !is.na(candidate_zip))

log_msg("Attaching ACS-derived ZCTA SVI proxy")
svi <- make_complete_acs_svi_proxy(community_path)
analysis_dat <- analysis_dat %>%
  left_join(svi, by = c("candidate_zip" = "zip", "listing_year_int" = "analysis_year"), suffix = c("", "_community")) %>%
  mutate(zcta_svi_proxy = coalesce(zcta_svi_proxy, zcta_svi_proxy_community)) %>%
  select(-any_of("zcta_svi_proxy_community")) %>%
  as.data.table()

log_msg("Reading exposure tables")
exposure_tables <- read_exposure_tables()
setkey(exposure_tables$pm25, zip, year, month)
setkey(exposure_tables$o3, zip, year, month)
setkey(exposure_tables$no2_monthly, zip, year, month)
setkey(exposure_tables$no2_annual, zip, year)

score_sources <- make_time_updated_score_sources()

pollutant_specs <- tribble(
  ~pollutant, ~exposure_end_date, ~term, ~label,
  "pm25", daily_pollution_end_date, "pm25_interval_5ug", "Monthly PM2.5 per 5 ug/m3",
  "o3", daily_pollution_end_date, "o3_interval_10ppb", "Monthly O3 per 10 ppb",
  "no2", as.Date("2025-12-31"), "no2_interval_10ppb", "Monthly NO2 per 10 ppb"
)

all_results <- list()
all_summaries <- list()
for (i in seq_len(nrow(pollutant_specs))) {
  spec <- pollutant_specs[i, ]
  log_msg("Building monthly time-varying severity intervals for ", spec$pollutant)
  intervals <- make_month_intervals(analysis_dat, spec$exposure_end_date)
  intervals <- add_time_updated_scores(intervals, score_sources)
  intervals[, pollutant := spec$pollutant]

  if (spec$pollutant == "pm25") {
    intervals <- exposure_tables$pm25[intervals, on = c("zip" = "candidate_zip", "year", "month")]
    intervals[, pm25_interval_5ug := pm25_interval_ug_m3 / 5]
  } else if (spec$pollutant == "o3") {
    intervals <- exposure_tables$o3[intervals, on = c("zip" = "candidate_zip", "year", "month")]
    intervals[, o3_interval_10ppb := o3_interval_ppb / 10]
  } else if (spec$pollutant == "no2") {
    intervals <- exposure_tables$no2_monthly[intervals, on = c("zip" = "candidate_zip", "year", "month")]
    intervals <- exposure_tables$no2_annual[intervals, on = c("zip", "year")]
    intervals[, no2_interval_ppb_final := fifelse(!is.na(no2_interval_ppb), no2_interval_ppb, no2_annual_ppb)]
    intervals[, no2_interval_10ppb := no2_interval_ppb_final / 10]
  }

  all_summaries[[spec$pollutant]] <- intervals[
    ,
    .(
      intervals = .N,
      candidate_organ_episodes = uniqueN(waitlist_row_id),
      people = uniqueN(PERS_ID),
      complete_exposure_intervals = sum(!is.na(get(spec$term))),
      min_year = min(year, na.rm = TRUE),
      max_year = max(year, na.rm = TRUE)
    ),
    by = .(pollutant, WL_ORG)
  ]

  for (org in target_organs) {
    result <- fit_tv_severity_model(intervals, org, spec$pollutant, spec$term, spec$label)
    all_results[[length(all_results) + 1L]] <- result
  }

  rm(intervals)
  invisible(gc())
}

tv_results <- bind_rows(all_results)
write_csv(tv_results, file.path(out_dir, "timevarying_pollution_severity_results.csv"))
write_csv(
  tv_results %>%
    mutate(
      estimate_ci = sprintf("%.3f (%.3f to %.3f)", estimate, conf_low, conf_high),
      p_value_display = if_else(p_value < 0.001, "<.001", sprintf("%.3f", p_value))
    ) %>%
    select(model_type, organ_label, outcome_label, pollutant, exposure, candidate_organ_episodes, intervals, estimate_ci, p_value_display, adjustment_set),
  file.path(out_dir, "timevarying_pollution_severity_table.csv")
)
write_csv(bind_rows(all_summaries), file.path(out_dir, "timevarying_pollution_severity_interval_summary.csv"))

log_msg("Wrote time-varying pollution severity outputs with listing year to ", normalizePath(out_dir, winslash = "/"))
